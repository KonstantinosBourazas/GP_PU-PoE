## AUC, rank and ESS summaries over the 100 datasets. Each correlation compares
## two settings over the same 100 datasets.
mc100_safe_cor <- function(x, y, method = "pearson") {
  x <- as.numeric(x)
  y <- as.numeric(y)
  mc100_assert(
    length(x) == length(y) &&
      length(x) >= 3L &&
      all(is.finite(x)) &&
      all(is.finite(y)),
    "Invalid paired correlation inputs."
  )
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::cor(x, y, method = method))
}

mc100_q <- function(x, p) unname(stats::quantile(x, probs = p, type = 7))

mc100_robustness_row <- function(z, baseline, id, r) {
  d <- z$node_summary
  b <- baseline$node_summary
  p <- d$latent_probability_mean[match(1:100, d$node)]
  p0 <- b$latent_probability_mean[match(1:100, b$node)]
  dg <- z$diagnostics
  data.frame(
    fit_id = id,
    replicate = as.integer(r),
    spearman_with_default = mc100_safe_cor(p, p0, "spearman"),
    spearman_unlabelled_with_default = mc100_safe_cor(
      p[1:80],
      p0[1:80],
      "spearman"
    ),
    mean_absolute_probability_difference = mean(abs(p - p0)),
    max_absolute_probability_difference = max(abs(p - p0)),
    mean_HPD80_width = mean(d$latent_HPD80_width),
    overall_AUC = dg$overall_AUC_latent_mean,
    hidden_AUC = dg$hidden_AUC_latent_mean,
    overall_LogLoss = dg$overall_LogLoss_latent_mean,
    overall_Brier = dg$overall_Brier_latent_mean,
    runtime_sec = dg$runtime_sec,
    stringsAsFactors = FALSE
  )
}

mc100_auc_reports <- function(robustness) {
  ids <- mc100_ids()
  reps <- 1:100
  mc100_assert(
    nrow(robustness) == 1300L &&
      !anyDuplicated(paste(robustness$fit_id, robustness$replicate)),
    "Expected 1,300 distinct setting/replication metric rows."
  )
  wide_for <- function(metric) {
    mat <- vapply(
      ids,
      function(id) {
        d <- robustness[robustness$fit_id == id, , drop = FALSE]
        mc100_assert(
          nrow(d) == 100L && setequal(d$replicate, reps),
          "Incomplete paired AUC data."
        )
        as.numeric(d[[metric]][match(reps, d$replicate)])
      },
      numeric(100)
    )
    mc100_assert(
      all(is.finite(mat)) && all(mat >= 0 & mat <= 1),
      "Invalid AUC values."
    )
    colnames(mat) <- ids
    mat
  }
  summary_rows <- pairs <- spreads <- spread_summaries <- matrices <- list()
  for (metric in c("overall_AUC", "hidden_AUC")) {
    mat <- wide_for(metric)
    baseline <- mat[, "DEFAULT"]
    for (id in ids) {
      x <- mat[, id]
      delta <- x - baseline
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        Metric = metric,
        Setting = id,
        Replications = 100L,
        Mean_AUC = mean(x),
        SD_AUC = stats::sd(x),
        Min_AUC = min(x),
        Q05_AUC = mc100_q(x, .05),
        Median_AUC = stats::median(x),
        Q95_AUC = mc100_q(x, .95),
        Max_AUC = max(x),
        Mean_difference_vs_DEFAULT = mean(delta),
        SD_difference_vs_DEFAULT = stats::sd(delta),
        MCSE_mean_difference = stats::sd(delta) / 10,
        Min_difference = min(delta),
        Q05_difference = mc100_q(delta, .05),
        Q95_difference = mc100_q(delta, .95),
        Max_difference = max(delta),
        Mean_absolute_difference = mean(abs(delta)),
        Max_absolute_difference = max(abs(delta)),
        Pearson_AUC_vs_DEFAULT = mc100_safe_cor(x, baseline),
        Spearman_AUC_vs_DEFAULT = mc100_safe_cor(x, baseline, "spearman"),
        stringsAsFactors = FALSE
      )
    }
    for (a in seq_len(length(ids) - 1L)) {
      for (b in seq.int(a + 1L, length(ids))) {
        pairs[[length(pairs) + 1L]] <- data.frame(
          Metric = metric,
          Setting_1 = ids[a],
          Setting_2 = ids[b],
          Replications = 100L,
          Pearson = mc100_safe_cor(mat[, a], mat[, b]),
          Spearman = mc100_safe_cor(mat[, a], mat[, b], "spearman"),
          stringsAsFactors = FALSE
        )
      }
    }
    for (method in c("pearson", "spearman")) {
      cm <- matrix(
        NA_real_,
        length(ids),
        length(ids),
        dimnames = list(ids, ids)
      )
      for (a in seq_along(ids)) {
        for (b in seq_along(ids)) {
          cm[a, b] <- mc100_safe_cor(mat[, a], mat[, b], method)
        }
      }
      matrices[[paste0(metric, "_", method)]] <- data.frame(
        Setting = ids,
        cm,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
    lo <- apply(mat, 1L, min)
    hi <- apply(mat, 1L, max)
    span <- hi - lo
    spreads[[metric]] <- data.frame(
      Metric = metric,
      replicate = reps,
      Min_AUC_across_settings = lo,
      Max_AUC_across_settings = hi,
      AUC_range_across_settings = span,
      stringsAsFactors = FALSE
    )
    means <- colMeans(mat)
    spread_summaries[[metric]] <- data.frame(
      Metric = metric,
      Replications = 100L,
      Settings = 13L,
      Minimum_setting_mean_AUC = min(means),
      Maximum_setting_mean_AUC = max(means),
      Range_of_setting_means = max(means) - min(means),
      Mean_within_dataset_AUC_range = mean(span),
      SD_within_dataset_AUC_range = stats::sd(span),
      Min_within_dataset_AUC_range = min(span),
      Q05_within_dataset_AUC_range = mc100_q(span, .05),
      Median_within_dataset_AUC_range = stats::median(span),
      Q95_within_dataset_AUC_range = mc100_q(span, .95),
      Max_within_dataset_AUC_range = max(span),
      stringsAsFactors = FALSE
    )
  }
  ranking <- list()
  for (id in ids) {
    for (target in c("All units", "Unlabelled units")) {
      d <- robustness[robustness$fit_id == id, , drop = FALSE]
      col <- if (target == "All units") {
        "spearman_with_default"
      } else {
        "spearman_unlabelled_with_default"
      }
      x <- d[[col]]
      mc100_assert(
        length(x) == 100L && all(is.finite(x)),
        "Missing within-dataset rank correlations."
      )
      ranking[[length(ranking) + 1L]] <- data.frame(
        Setting = id,
        Target = target,
        Replications = 100L,
        Mean_Spearman = mean(x),
        SD_Spearman = stats::sd(x),
        Min_Spearman = min(x),
        Q05_Spearman = mc100_q(x, .05),
        Median_Spearman = stats::median(x),
        Q95_Spearman = mc100_q(x, .95),
        Max_Spearman = max(x),
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    summary = mc100_bind(summary_rows),
    pairwise_correlations = mc100_bind(pairs),
    correlation_matrices = matrices,
    ranges_by_replication = mc100_bind(spreads),
    range_summary = mc100_bind(spread_summaries),
    ranking_stability = mc100_bind(ranking)
  )
}

mc100_save_auc_reports <- function(reports, folder) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  tab <- function(x, name) {
    utils::write.csv(x, file.path(folder, name), row.names = FALSE, na = "NA")
  }
  tab(reports$summary, "MC100_AUC_summary_and_paired_differences.csv")
  tab(reports$pairwise_correlations, "MC100_AUC_pairwise_correlations.csv")
  for (name in names(reports$correlation_matrices)) {
    tab(
      reports$correlation_matrices[[name]],
      paste0("MC100_AUC_correlation_matrix_", name, ".csv")
    )
  }
  tab(reports$ranges_by_replication, "MC100_AUC_ranges_by_replication.csv")
  tab(reports$range_summary, "MC100_AUC_range_summary.csv")
  tab(reports$ranking_stability, "MC100_ranking_stability_vs_DEFAULT.csv")
  invisible(TRUE)
}

mc100_original_node_variation <- function(results, id) {
  mc100_assert(
    length(results) == 100L,
    "Expected 100 datasets for node variation."
  )
  mat <- vapply(
    results,
    function(z) {
      d <- z$node_summary
      as.numeric(d$latent_probability_mean[match(1:100, d$node)])
    },
    numeric(100)
  )
  mc100_assert(
    all(is.finite(mat)),
    "Missing posterior means for node variation."
  )
  data.frame(
    fit_id = id,
    original_node = 1:100,
    Replications = 100L,
    Mean_posterior_mean = rowMeans(mat),
    SD_posterior_mean_across_datasets = apply(mat, 1L, stats::sd),
    MCSE_average_posterior_mean = apply(mat, 1L, stats::sd) / 10,
    Q05_posterior_mean = apply(mat, 1L, mc100_q, p = .05),
    Q95_posterior_mean = apply(mat, 1L, mc100_q, p = .95),
    stringsAsFactors = FALSE
  )
}

## ESS of the fitted probabilities over the 100 datasets. Reads only the ESS
## summaries saved in the 1,300 checkpoints; no ESS is recomputed from chains.
mc100_ess_summaries <- function(archive) {
  columns <- c(
    "min_ESS",
    "q25_ESS",
    "median_ESS",
    "mean_ESS",
    "q75_ESS",
    "max_ESS"
  )
  default_hash <- vapply(
    1:100,
    function(r) {
      z <- mc100_read(mc100_path(archive, "DEFAULT", r))
      mc100_assert(
        is.list(z) &&
          isTRUE(z$complete) &&
          identical(z$fit_id, "DEFAULT") &&
          identical(as.integer(z$replicate), as.integer(r)) &&
          length(z$toy_hash) == 1L,
        paste0("Invalid DEFAULT checkpoint for replication ", r, ".")
      )
      as.character(z$toy_hash)
    },
    character(1)
  )
  rows <- vector("list", 1300L)
  k <- 0L
  for (id in mc100_ids()) {
    for (r in 1:100) {
      path <- mc100_path(archive, id, r)
      z <- mc100_read(path)
      tag <- paste0(id, " replication ", r)
      mc100_assert(
        is.list(z) &&
          isTRUE(z$complete) &&
          identical(z$fit_id, id) &&
          identical(as.integer(z$replicate), as.integer(r)),
        paste0(tag, ": checkpoint identity or completion mismatch.")
      )
      mc100_assert(
        identical(as.character(z$toy_hash), default_hash[r]),
        paste0(tag, ": the dataset differs from DEFAULT.")
      )
      es <- z$node_ess_summary
      mc100_assert(
        is.data.frame(es) && nrow(es) == 1L && all(columns %in% names(es)),
        paste0(tag, ": invalid node_ess_summary.")
      )
      values <- as.numeric(es[1L, columns])
      mc100_assert(
        all(is.finite(values)) && all(values > 0),
        paste0(tag, ": invalid ESS values.")
      )
      ne <- z$node_probability_ess
      mc100_assert(
        is.data.frame(ne) &&
          nrow(ne) == 100L &&
          all(c("ESS", "retained_draws") %in% names(ne)) &&
          all(is.finite(ne$ESS)) &&
          all(ne$ESS > 0) &&
          all(ne$retained_draws == 5000L),
        paste0(tag, ": invalid node_probability_ess.")
      )
      mc100_assert(
        isTRUE(all.equal(
          stats::median(ne$ESS),
          as.numeric(es$median_ESS),
          tolerance = 1e-8
        )) &&
          isTRUE(all.equal(
            mean(ne$ESS),
            as.numeric(es$mean_ESS),
            tolerance = 1e-8
          )),
        paste0(tag, ": the ESS summary does not match the 100 node ESS values.")
      )
      k <- k + 1L
      rows[[k]] <- data.frame(
        Setting = id,
        Replication = r,
        Median_ESS = as.numeric(es$median_ESS),
        Mean_ESS = as.numeric(es$mean_ESS),
        Minimum_ESS = as.numeric(es$min_ESS),
        Q25_ESS = as.numeric(es$q25_ESS),
        Q75_ESS = as.numeric(es$q75_ESS),
        Maximum_ESS = as.numeric(es$max_ESS),
        stringsAsFactors = FALSE
      )
    }
  }
  by_replication <- mc100_bind(rows)
  summary <- mc100_bind(lapply(mc100_ids(), function(id) {
    d <- by_replication[by_replication$Setting == id, , drop = FALSE]
    mc100_assert(
      nrow(d) == 100L,
      paste0("Expected 100 replications for ", id, ".")
    )
    data.frame(
      Setting = id,
      Replications = nrow(d),
      Mean_of_median_ESS = mean(d$Median_ESS),
      SD_of_median_ESS = stats::sd(d$Median_ESS),
      Mean_of_mean_ESS = mean(d$Mean_ESS),
      SD_of_mean_ESS = stats::sd(d$Mean_ESS),
      Median_of_median_ESS = stats::median(d$Median_ESS),
      Median_of_mean_ESS = stats::median(d$Mean_ESS),
      Minimum_median_ESS = min(d$Median_ESS),
      Minimum_mean_ESS = min(d$Mean_ESS),
      Minimum_node_ESS = min(d$Minimum_ESS),
      stringsAsFactors = FALSE
    )
  }))
  list(by_replication = by_replication, summary = summary)
}
