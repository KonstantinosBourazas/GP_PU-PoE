## Checks of the archived results.
sar_unpack_methods <- function(z) {
  if (!is.null(z$method_results)) unname(z$method_results) else list(z)
}

sar_validate_method <- function(v, r, group, expected_y = NULL) {
  m <- v$metrics
  d <- v$node_scores
  tag <- paste(group, r, if (is.data.frame(m)) m$method[1L] else "unknown")
  check <- function(ok, text) sar_assert(ok, paste(tag, text))
  check(
    is.data.frame(m) && nrow(m) == 1L && is.data.frame(d) && nrow(d) == 400L,
    "invalid result dimensions"
  )
  check(
    all(m$replicate == r) &&
      all(d$replicate == r) &&
      !anyDuplicated(d$node) &&
      setequal(d$node, 1:400),
    "replication/node mismatch"
  )
  d <- d[match(1:400, d$node), , drop = FALSE]
  truth <- c(rep(0L, 355), rep(1L, 45))
  observed <- as.integer(d$Y)
  hidden <- which(truth == 1L & observed == 0L)
  unlabelled <- which(observed == 0L)
  check(
    identical(as.integer(d$T), truth) &&
      all(observed %in% c(0L, 1L)) &&
      all(observed <= truth) &&
      sum(observed) == 30L &&
      length(hidden) == 15L,
    "invalid SAR labels: the hidden indices must be determined separately in each dataset"
  )
  if (!is.null(expected_y)) {
    check(
      identical(observed, as.integer(expected_y)),
      "observed labels differ from the common SAR dataset"
    )
  }
  if ("is_hidden_positive" %in% names(d)) {
    check(
      identical(
        as.integer(d$is_hidden_positive),
        as.integer(seq_len(400L) %in% hidden)
      ),
      "hidden-positive indicator disagrees with T and Y"
    )
  }
  p <- as.numeric(d$probability)
  score <- if ("decision_score" %in% names(d)) {
    as.numeric(d$decision_score)
  } else {
    p
  }
  check(all(is.finite(c(p, score))) && all(p >= 0 & p <= 1), "invalid scores")
  rk <- rank(-score, ties.method = "average")
  ru <- rank(-score[unlabelled], ties.method = "average")
  clipped <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  point <- c(
    Overall_AUC = sar_auc(score, truth),
    Hidden_AUC = sar_auc(score[unlabelled], truth[unlabelled]),
    Overall_Recall_at_45 = sum(truth == 1L & rk <= 45) / 45,
    Hidden_Recall_at_15 = sum(truth[unlabelled] == 1L & ru <= 15) / 15,
    Overall_LogLoss = -mean(
      truth * log(clipped) + (1 - truth) * log(1 - clipped)
    ),
    Hidden_LogLoss = -mean(log(clipped[hidden])),
    Overall_Brier = mean((p - truth)^2),
    Hidden_Brier = mean((p[hidden] - 1)^2)
  )
  check(
    all(names(point) %in% names(m)) &&
      sar_near(point, unlist(m[1L, names(point)]), 1e-8),
    "point metrics disagree with the saved node scores"
  )
  for (level in c(80L, 95L)) {
    lower <- paste0("PI", level, "_lower")
    upper <- paste0("PI", level, "_upper")
    cols <- paste0(
      c(
        "Overall_Coverage_",
        "Hidden_Coverage_",
        "Overall_Interval_Score_",
        "Hidden_Interval_Score_"
      ),
      level
    )
    check(
      all(c(lower, upper) %in% names(d)) && all(cols %in% names(m)),
      "missing predictive-interval fields"
    )
    lo <- d[[lower]]
    hi <- d[[upper]]
    if (all(is.na(lo)) && all(is.na(hi))) {
      check(
        all(is.na(unlist(m[1L, cols]))),
        "unavailable uncertainty must stay NA"
      )
      next
    }
    check(
      !anyNA(c(lo, hi)) &&
        all(lo %in% c(0, 1)) &&
        all(hi %in% c(0, 1)) &&
        all(lo <= hi),
      "invalid binary predictive interval"
    )
    covered <- as.numeric(lo <= truth & truth <= hi)
    alpha <- 1 - level / 100
    score_is <- hi -
      lo +
      2 / alpha * (lo - truth) * (truth < lo) +
      2 / alpha * (truth - hi) * (truth > hi)
    calculated <- c(
      mean(covered),
      mean(covered[hidden]),
      mean(score_is),
      mean(score_is[hidden])
    )
    check(
      sar_near(calculated, unlist(m[1L, cols]), 1e-8),
      "coverage/interval scores disagree with saved node limits"
    )
  }
  m$Method <- sar_canonical(m$method)
  m$group <- group
  d$Method <- sar_canonical(d$method)
  d$group <- group
  list(metrics = m, scores = d)
}

sar_validate_selection <- function(x, r, expected_y) {
  s <- x$sar_selection_summary
  d <- x$sar_selection_diagnostics
  tag <- paste("SAR selection, replication", r)
  sar_assert(
    is.data.frame(s) && nrow(s) == 1L && is.data.frame(d) && nrow(d) == 45L,
    paste(tag, "missing selection diagnostics")
  )
  sar_assert(
    s$lambda_SAR == 2 &&
      s$z_clip == 3 &&
      s$n_hidden == 15L &&
      s$n_observed_positive == 30L &&
      s$sar_seed == 20260718 + 350000 + 1000 * r + 71,
    paste(tag, "wrong lambda/counts/seed")
  )
  sar_assert(
    !anyDuplicated(d$node) && setequal(d$node, 356:400),
    paste(tag, "invalid positive indices")
  )
  d <- d[match(356:400, d$node), , drop = FALSE]
  hidden <- which(expected_y == 0L & seq_along(expected_y) >= 356L)
  sar_assert(
    identical(as.integer(d$hidden_SAR), as.integer(d$node %in% hidden)) &&
      identical(as.integer(d$observed_SAR), as.integer(!d$node %in% hidden)),
    paste(tag, "indicators do not match labels")
  )
  z <- as.numeric(d$standardized_oracle_score)
  w <- exp(-2 * z - max(-2 * z))
  sar_assert(
    all(is.finite(z)) &&
      all(abs(z) <= 3 + 1e-12) &&
      sar_near(w, d$hidden_weight, 1e-10),
    paste(tag, "hiding weights do not follow the saved SAR rule")
  )
  h <- d$hidden_SAR == 1L
  scores <- d$oracle_covariate_score
  sar_assert(
    sar_near(
      c(mean(scores[h]), mean(scores[!h]), mean(scores[!h]) - mean(scores[h])),
      unlist(s[
        1L,
        c(
          "mean_score_hidden",
          "mean_score_observed",
          "observed_minus_hidden_score"
        )
      ]),
      1e-9
    ),
    paste(tag, "selection summaries disagree with saved scores")
  )
  s$replicate <- r
  s
}

sar_read_group <- function(archive, group, expected_labels, reps = 1:100) {
  dir <- file.path(archive, group)
  sp <- sar_group_spec()
  sp <- sp[sp$group == group, , drop = FALSE]
  paths <- file.path(dir, sp$checkpoint_dir, sprintf(sp$pattern, reps))
  expected <- sar_registry()
  expected <- expected$Method[expected$group == group]
  out <- vector("list", length(reps))
  selections <- vector("list", length(reps))
  sig <- NULL
  for (j in seq_along(reps)) {
    r <- reps[j]
    x <- sar_read(paths[j])
    sar_assert(
      isTRUE(x$complete) && identical(as.integer(x$replicate), as.integer(r)),
      paste(group, r, "incomplete checkpoint")
    )
    if (is.null(sig)) {
      sig <- x$settings_signature
    }
    sar_assert(
      identical(x$settings_signature, sig),
      paste(group, r, "different settings")
    )
    lm <- sig$label_mechanism
    sar_assert(
      is.list(lm) &&
        lm$lambda == 2 &&
        lm$z_clip == 3 &&
        lm$exactly_n_hidden_without_replacement == 15L,
      paste(group, r, "not the SAR design with lambda = 2")
    )
    geo <- x$geometry_summary
    expected_seeds <- c(
      covariate_seed = 20260718 + 100000 + 1000 * r + 11,
      network_seed = 20260718 + 200000 + 1000 * r + 29,
      modularity_seed = 20260718 + 300000 + 1000 * r + 43,
      sar_seed = 20260718 + 350000 + 1000 * r + 71
    )
    present <- intersect(names(expected_seeds), names(geo))
    sar_assert(
      all(c("covariate_seed", "sar_seed") %in% present) &&
        sar_near(
          unlist(geo[1L, present, drop = FALSE]),
          expected_seeds[present],
          0
        ),
      paste(group, r, "different seed plan")
    )
    if (!group %in% c("Competitors", "nnPU")) {
      sar_assert(
        sar_near(
          unlist(sig[c(
            "production_burnin",
            "production_sampling",
            "production_thin"
          )]),
          c(15000, 40000, 20),
          0
        ),
        paste(group, "not the full schedule")
      )
      saved <- vapply(
        sar_unpack_methods(x),
        function(v) as.numeric(v$n_saved),
        numeric(1)
      )
      sar_assert(all(saved == 2000L), paste(group, r, "retained draw count"))
    }
    if (group == "Bayesian_linear") {
      ev <- sig$evaluation$predictive_label_draws
      sar_assert(
        is.character(ev) &&
          length(ev) == 1L &&
          grepl("Bernoulli(plogis", ev, fixed = TRUE),
        "Unexpected evaluation version of the Bayesian linear models."
      )
      ap <- sig$priors
      lp <- sar_linear_priors()
      sar_validate_linear_priors(ap, "Bayesian linear archive")
    }
    pp <- lapply(
      sar_unpack_methods(x),
      sar_validate_method,
      r = r,
      group = group,
      expected_y = expected_labels[[r]]
    )
    got <- unlist(lapply(pp, function(v) v$metrics$Method), use.names = FALSE)
    sar_assert(
      !anyDuplicated(got) && setequal(got, expected),
      paste(group, r, "missing/extra methods")
    )
    out[[j]] <- pp
    selections[[j]] <- sar_validate_selection(x, r, expected_labels[[r]])
  }
  flat <- do.call(c, out)
  m <- sar_bind(lapply(flat, `[[`, "metrics"))
  d <- sar_bind(lapply(flat, `[[`, "scores"))
  csv <- utils::read.csv(
    file.path(dir, paste0(sp$prefix, "_metrics_by_replicate.csv")),
    stringsAsFactors = FALSE
  )
  for (method in unique(m$method)) {
    aa <- m[m$method == method, , drop = FALSE]
    bb <- csv[csv$method == method, , drop = FALSE]
    aa <- aa[order(aa$replicate), , drop = FALSE]
    bb <- bb[order(bb$replicate), , drop = FALSE]
    sar_assert(
      nrow(aa) == 100L &&
        nrow(bb) == 100L &&
        !anyDuplicated(bb$replicate) &&
        setequal(bb$replicate, 1:100),
      paste(group, method, "CSV count")
    )
    for (metric in sar_metrics()) {
      sar_assert(
        sar_near(aa[[metric]], bb[[metric]], 1e-8),
        paste(group, method, metric, "CSV/checkpoint mismatch")
      )
    }
  }
  selection <- sar_bind(selections)
  selection$group <- group
  cat(
    "PASS: ",
    group,
    ": 100 replications, ",
    length(expected),
    if (length(expected) == 1L) " method" else " methods",
    ", shared SAR labels.\n",
    sep = ""
  )
  list(metrics = m, scores = d, selection = selection, signature = sig)
}

sar_archive_manifest_check <- function(root, archive) {
  f <- utils::read.csv(
    file.path(sar_source_root(root), "reference/ARCHIVE_MANIFEST.csv"),
    stringsAsFactors = FALSE
  )
  paths <- file.path(archive, f$file)
  missing <- !file.exists(paths)
  sar_assert(
    !any(missing),
    paste(
      "Missing archive files:",
      paste(head(f$file[missing], 8L), collapse = "\n")
    )
  )
  sar_assert(
    identical(as.character(sar_md5(paths)), as.character(f$md5)),
    "Archive files differ from ARCHIVE_MANIFEST.csv."
  )
  cat("PASS: ", nrow(f), " archive files match the manifest.\n", sep = "")
  inventory <- utils::read.csv(
    file.path(
      root,
      "simulation_nonlinear_SAR_all18/reference/ARCHIVE_MANIFEST.csv"
    ),
    stringsAsFactors = FALSE
  )
  expected_linear <- sub(
    "^Bayesian_linear/",
    "",
    inventory$file[startsWith(inventory$file, "Bayesian_linear/")]
  )
  actual_linear <- list.files(
    file.path(archive, "Bayesian_linear"),
    recursive = TRUE,
    all.files = TRUE
  )
  sar_assert(
    setequal(actual_linear, expected_linear) &&
      length(actual_linear) == length(expected_linear),
    "Unexpected extra or missing Bayesian archive files"
  )
  invisible(TRUE)
}

sar_collect_archive <- function(root, archive) {
  sar_linear_archive_checks(archive)

  sar_archive_manifest_check(root, archive)
  sar_pilot_archive_checks(archive)
  labels <- vector("list", 100L)
  for (r in 1:100) {
    fn <- sprintf("replicate_%03d_geometry.rds", r)
    a <- sar_read(file.path(archive, "Competitors/generated_geometries", fn))
    b <- sar_read(file.path(archive, "nnPU/generated_geometries", fn))
    for (field in c(
      "X_raw",
      "X_std",
      "A_network",
      "Z_network_modularity",
      "T",
      "Y",
      "hidden_positive_idx",
      "observed_positive_idx"
    )) {
      sar_assert(
        !is.null(a[[field]]) &&
          !is.null(b[[field]]) &&
          isTRUE(all.equal(
            a[[field]],
            b[[field]],
            tolerance = 0,
            check.attributes = FALSE
          )),
        paste("Different saved competitor/nnPU dataset", r, field)
      )
    }
    labels[[r]] <- as.integer(a$Y)
  }
  cat(
    "PASS: 100 datasets shared by the competitors and nnPU; the other methods have the same Y, T and seeds.\n"
  )
  groups <- lapply(sar_groups(), function(g) sar_read_group(archive, g, labels))
  names(groups) <- sar_groups()
  metrics <- sar_bind(lapply(groups, `[[`, "metrics"))
  scores <- sar_bind(lapply(groups, `[[`, "scores"))
  selection <- sar_bind(lapply(groups, `[[`, "selection"))
  sar_assert(
    nrow(metrics) == 1800L && nrow(scores) == 720000L,
    "Incomplete all-method archive."
  )
  list(
    metrics = metrics,
    scores = scores,
    selection = selection,
    mode = "archived_full",
    replications = 100L
  )
}

sar_pilot_archive_checks <- function(archive) {
  specs <- sar_group_spec()
  total <- 0L
  for (group in setdiff(sar_groups(), c("Competitors", "nnPU"))) {
    folder <- file.path(archive, group)
    files <- list.files(
      file.path(folder, "pilot_checkpoints"),
      pattern = "[.]rds$",
      full.names = TRUE
    )
    expected <- if (group == "Bayesian_linear") 20L else 5L
    sar_assert(
      length(files) == expected,
      paste(group, "pilot checkpoint count")
    )
    for (path in files) {
      x <- sar_read(path)
      sig <- x$settings_signature
      sar_assert(
        is.list(x$tail) &&
          is.list(x$final_state) &&
          (group == "Bayesian_linear" || isTRUE(x$complete)),
        paste("Incomplete pilot", path)
      )
      sar_assert(
        sig$pilot_burnin == 40000L &&
          sig$pilot_tail_length == 5000L &&
          sig$label_mechanism$lambda == 2 &&
          sig$label_mechanism$z_clip == 3,
        paste("Pilot settings mismatch", path)
      )
      if (group == "Bayesian_linear") {
        ap <- sig$priors
        sar_validate_linear_priors(ap, "Bayesian linear archive")
      }
    }
    sp <- specs[specs$group == group, , drop = FALSE]
    stem <- if (group == "Bayesian_linear") {
      "global_pilot_calibration_5_datasets"
    } else {
      paste0(sp$prefix, "_global_pilot_calibration_5_datasets")
    }
    hash <- sar_md5(file.path(folder, paste0(stem, ".rds")))
    first <- sar_read(file.path(
      folder,
      sp$checkpoint_dir,
      sprintf(sp$pattern, 1L)
    ))
    sar_assert(
      identical(hash, first$settings_signature$global_calibration_hash),
      paste(group, "production does not reference this pilot calibration")
    )
    total <- total + length(files)
  }
  cat(
    "PASS: ",
    total,
    " pilot checkpoints and six referenced global calibrations.\n",
    sep = ""
  )
  invisible(TRUE)
}
