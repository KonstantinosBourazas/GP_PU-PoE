lscar_summarize <- function(metrics, reps) {
  r <- lscar_registry()
  out <- list()
  i <- 0L
  for (method in r$Method) {
    d <- metrics[metrics$Method == method, , drop = FALSE]
    lscar_assert(
      nrow(d) == length(reps) &&
        !anyDuplicated(d$replicate) &&
        setequal(d$replicate, reps),
      paste("Incomplete replications:", method)
    )
    for (metric in lscar_metrics()) {
      x <- as.numeric(d[[metric]])
      available <- !all(is.na(x))
      lscar_assert(
        !available || all(is.finite(x)),
        paste("Partial missing metric:", method, metric)
      )
      i <- i + 1L
      out[[i]] <- data.frame(
        Method = method,
        Metric = metric,
        n = length(reps),
        mean = if (available) mean(x) else NA_real_,
        sd = if (available) stats::sd(x) else NA_real_,
        mcse = if (available) stats::sd(x) / sqrt(length(x)) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  ans <- do.call(rbind, out)
  ans$formatted <- lscar_format(ans$mean, ans$sd)
  ans
}

lscar_wide_table <- function(summary) {
  methods <- lscar_registry()$Method
  out <- data.frame(Method = methods, stringsAsFactors = FALSE)
  for (metric in lscar_metrics()) {
    d <- summary[summary$Metric == metric, ]
    out[[metric]] <- d$formatted[match(methods, d$Method)]
  }
  out
}

lscar_eta_table <- function(metrics) {
  ids <- c("GP-PU-PoE", "GP-PU-PoE_beta0")
  lscar_bind(lapply(ids, function(id) {
    d <- metrics[metrics$Method == id, , drop = FALSE]
    p <- d$eta_posterior_mean
    lo <- d$eta_q10
    hi <- d$eta_q90
    truth <- 1 / 3
    lscar_assert(
      length(p) > 0L && all(is.finite(c(p, lo, hi))) && all(lo <= hi),
      paste("Missing eta quantiles:", id)
    )
    H <- 30 * p / (1 - p)
    data.frame(
      Method = id,
      Replications = nrow(d),
      Mean_eta_hat = mean(p),
      RMSE_eta_hat = sqrt(mean((p - truth)^2)),
      Coverage_80_eta = mean(lo <= truth & truth <= hi),
      Mean_IS_80_eta = mean(
        hi -
          lo +
          10 * (lo - truth) * (truth < lo) +
          10 * (truth - hi) * (truth > hi)
      ),
      Mean_H_hat_eta = mean(H),
      RMSE_H_hat_eta = sqrt(mean((H - 15)^2)),
      stringsAsFactors = FALSE
    )
  }))
}

## Compares the tables with the reference values in reference/EXPECTED_*.csv.
lscar_check_reference <- function(root, summary, eta) {
  near <- function(a, b, tol = 1e-9) {
    length(a) == length(b) &&
      all(is.na(a) == is.na(b)) &&
      all(abs(a[!is.na(a)] - b[!is.na(b)]) <= tol)
  }
  fmt <- function(m, s) {
    ifelse(is.na(m) | is.na(s), "--", sprintf("%.3f (%.3f)", m, s))
  }
  rd <- file.path(root, "simulation_linear_SCAR_all18", "reference")
  ref <- utils::read.csv(
    file.path(rd, "EXPECTED_ARCHIVED_SUMMARIES.csv"),
    stringsAsFactors = FALSE
  )
  key <- function(d) paste(d$Method, d$Metric, sep = "::")
  lscar_assert(
    nrow(ref) == 216L &&
      nrow(summary) == 216L &&
      !anyDuplicated(key(ref)) &&
      !anyDuplicated(key(summary)),
    "Invalid reference table"
  )
  k <- match(key(ref), key(summary))
  lscar_assert(!anyNA(k), "Missing cells in the comparison")
  ok <- vapply(
    seq_len(nrow(ref)),
    function(i) {
      near(
        c(summary$mean[k[i]], summary$sd[k[i]]),
        c(ref$mean[i], ref$sd[i]),
        1e-9
      )
    },
    logical(1)
  )
  out <- data.frame(
    Method = ref$Method,
    Metric = ref$Metric,
    computed = fmt(summary$mean[k], summary$sd[k]),
    reference = fmt(ref$mean, ref$sd),
    numeric_check_pass = ok,
    stringsAsFactors = FALSE
  )
  er <- utils::read.csv(
    file.path(rd, "EXPECTED_ETA_SUMMARIES.csv"),
    stringsAsFactors = FALSE
  )
  em <- c(
    "Mean_eta_hat",
    "RMSE_eta_hat",
    "Coverage_80_eta",
    "Mean_IS_80_eta",
    "Mean_H_hat_eta",
    "RMSE_H_hat_eta"
  )
  lscar_assert(
    nrow(er) == 2L &&
      !anyDuplicated(er$Method) &&
      setequal(er$Method, eta$Method),
    "Invalid reference table for eta"
  )
  for (i in seq_len(nrow(er))) {
    for (metric in em) {
      value <- eta[match(er$Method[i], eta$Method), metric]
      want <- er[i, metric]
      out <- rbind(
        out,
        data.frame(
          Method = er$Method[i],
          Metric = metric,
          computed = sprintf("%.3f", value),
          reference = sprintf("%.3f", want),
          numeric_check_pass = near(value, want, 1e-9),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  out$reference_match <- out$computed == out$reference
  out$pass <- out$numeric_check_pass
  out$status <- ifelse(
    !out$pass,
    "NUMERIC_MISMATCH",
    ifelse(out$reference_match, "MATCH", "MATCH_WITHIN_NUMERIC_PRECISION")
  )
  out
}

lscar_latex_table <- function(summary, eta, quick = FALSE) {
  rn <- lscar_registry()$Method
  lab <- function(x) ifelse(x == "GP-PU-PoE_beta0", "GP-PU-PoE$_{\\beta_0}$", x)
  metrics <- lscar_metrics()
  titles <- c(
    "Overall\\\\AUC",
    "Hidden\\\\AUC",
    "Overall\\\\Recall@45",
    "Hidden\\\\Recall@15",
    "Overall\\\\LogLoss",
    "Hidden\\\\LogLoss",
    "Overall\\\\Cov$_{80}$",
    "Hidden\\\\Cov$_{80}$",
    "Overall\\\\IS$_{80}$",
    "Hidden\\\\IS$_{80}$",
    "Overall\\\\Brier",
    "Hidden\\\\Brier"
  )
  best <- setNames(vector("list", length(metrics)), metrics)
  for (m in metrics) {
    d <- summary[summary$Metric == m & is.finite(summary$mean), ]
    objective <- round(d$mean, 3)
    if (grepl("Coverage", m)) {
      objective <- abs(objective - .80)
    } else if (grepl("AUC|Recall", m)) {
      objective <- -objective
    }
    best[[m]] <- d$Method[objective == min(objective)]
  }
  cell <- function(method, metric) {
    d <- summary[summary$Method == method & summary$Metric == metric, ]
    if (is.na(d$mean)) {
      return("--")
    }
    mn <- sprintf("%.3f", d$mean)
    if (!quick && method %in% best[[metric]]) {
      mn <- paste0("\\textbf{", mn, "}")
    }
    paste0(mn, " (", sprintf("%.3f", d$sd), ")")
  }
  lines <- c(
    "% Tables of the linear SCAR study, computed from the results of each replication.",
    "% Requires the LaTeX packages graphicx, booktabs and arydshln.",
    "\\begin{table}[!t]",
    "\\centering",
    paste0(
      "\\caption{",
      if (quick) "Quick check, not the results of the paper. " else "",
      "Linear design under SCAR. Monte Carlo means over ",
      unique(summary$n),
      " replications, standard deviations in parentheses.",
      if (!quick) {
        " Best value in bold, closest to nominal for coverage."
      } else {
        ""
      },
      "}"
    ),
    paste0("\\label{tab:sim-linear-scar", if (quick) "-quick" else "", "}"),
    "\\begingroup",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{1.4pt}",
    "\\renewcommand{\\arraystretch}{1.1}",
    "\\setlength{\\aboverulesep}{0.15ex}",
    "\\setlength{\\belowrulesep}{0.15ex}"
  )
  for (block in 1:2) {
    ix <- if (block == 1) 1:6 else 7:12
    lines <- c(
      lines,
      "\\resizebox{\\textwidth}{!}{",
      "\\begin{tabular}{lcccccc}",
      "\\toprule",
      paste0(
        "Method & ",
        paste(paste0("\\shortstack{", titles[ix], "}"), collapse = " & "),
        " ",
        "\\\\"
      )
    )
    lines <- c(lines, "\\midrule")
    for (j in seq_along(rn)) {
      lines <- c(
        lines,
        paste0(
          lab(rn[j]),
          " & ",
          paste(
            vapply(metrics[ix], function(m) cell(rn[j], m), character(1)),
            collapse = " & "
          ),
          " ",
          "\\\\"
        )
      )
      if (j %in% c(4, 9)) lines <- c(lines, "\\hdashline")
    }
    lines <- c(lines, "\\bottomrule", "\\end{tabular}", "}")
    if (block == 1) lines <- c(lines, "\\vspace{0.3ex}")
  }
  lines <- c(
    lines,
    "\\endgroup",
    "\\end{table}",
    "",
    "\\begin{table}[!t]",
    "\\centering",
    paste0(
      "\\caption{",
      if (quick) "Quick check, not the results of the paper. " else "",
      "Estimation of $\\eta$ and $H$ in the linear design under SCAR over ",
      eta$Replications[1L],
      " replications.}"
    ),
    paste0("\\label{tab:eta-linear-scar", if (quick) "-quick" else "", "}"),
    "\\small",
    "\\setlength{\\tabcolsep}{6pt}",
    "\\begin{tabular}{lrrrrrr}",
    "\\toprule",
    paste0(
      "Method & Mean $\\hat\\eta$ & RMSE$_\\eta$ & Cov$_{80,\\eta}$ & IS$_{80,\\eta}$ & Mean $\\hat H$ & RMSE$_H$ ",
      "\\\\"
    ),
    "\\midrule"
  )
  em <- c(
    "Mean_eta_hat",
    "RMSE_eta_hat",
    "Coverage_80_eta",
    "Mean_IS_80_eta",
    "Mean_H_hat_eta",
    "RMSE_H_hat_eta"
  )
  for (j in seq_len(nrow(eta))) {
    values <- vapply(
      em,
      function(m) {
        value <- sprintf("%.3f", eta[[m]][j])
        target <- switch(
          m,
          Mean_eta_hat = 1 / 3,
          Coverage_80_eta = .8,
          Mean_H_hat_eta = 15,
          0
        )
        winner <- which.min(abs(eta[[m]] - target))
        if (!quick && j == winner) paste0("\\textbf{", value, "}") else value
      },
      character(1)
    )
    lines <- c(
      lines,
      paste0(
        lab(eta$Method[j]),
        " & ",
        paste(values, collapse = " & "),
        " ",
        "\\\\"
      )
    )
  }
  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

lscar_write_reports <- function(root, bundle, out) {
  quick <- bundle$mode == "quick"
  suffix <- if (quick) "_QUICK" else ""
  summary <- lscar_summarize(bundle$metrics, if (quick) 1:2 else 1:100)
  eta <- lscar_eta_table(bundle$metrics)
  wide <- lscar_wide_table(summary)
  td <- file.path(out, "tables")
  vd <- file.path(out, "verification")
  lscar_csv(
    bundle$metrics,
    file.path(td, paste0("metrics_by_replication_all18", suffix, ".csv"))
  )
  lscar_csv(
    bundle$scores,
    file.path(td, paste0("node_scores_all18", suffix, ".csv"))
  )
  lscar_csv(
    summary,
    file.path(td, paste0("metric_summary_long", suffix, ".csv"))
  )
  lscar_csv(
    wide,
    file.path(td, paste0("Table_linear_SCAR_all18", suffix, ".csv"))
  )
  lscar_csv(eta, file.path(td, paste0("Table_eta_H", suffix, ".csv")))
  writeLines(
    lscar_latex_table(summary, eta, quick),
    file.path(td, paste0("Tables_linear_SCAR_all18", suffix, ".tex"))
  )
  if (bundle$mode == "archived_full") {
    comparison <- lscar_check_reference(root, summary, eta)
    lscar_csv(
      comparison,
      file.path(out, "verification", "comparison_results.csv")
    )
    lscar_assert(
      all(comparison$numeric_check_pass),
      "The comparison with the reference values failed. See verification/comparison_results.csv"
    )
    cat(
      "\nPASS: the 216 cells of the main table and the 12 eta and H values agree with the reference values.\n"
    )
  }
  print(wide, row.names = FALSE)
  cat("\nETA / H\n")
  print(eta, row.names = FALSE)
  lscar_atomic_rds(
    list(
      mode = bundle$mode,
      metrics = bundle$metrics,
      summary = summary,
      eta = eta
    ),
    file.path(out, "combined_results.rds")
  )
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  invisible(list(summary = summary, table = wide, eta = eta))
}
