## Checks the failure records of the archived checkpoints. Competitors,
## Bayesian_linear and nnPU store method_failures as a list, not as a data
## frame, and nrow() of an empty list is NULL. Real failure entries and
## unexpected types are rejected.
lscar_assert_empty_method_failures <- function(failures, tag) {
  if (is.null(failures)) {
    return(invisible(TRUE))
  }
  if (is.data.frame(failures)) {
    n_failures <- nrow(failures)
  } else if (
    is.list(failures) && !is.object(failures) && is.null(dim(failures))
  ) {
    n_failures <- length(failures)
  } else {
    stop(
      paste0(
        tag,
        " invalid method_failures format: expected NULL, a list ",
        "or a data frame."
      ),
      call. = FALSE
    )
  }
  lscar_assert(
    n_failures == 0L,
    paste0(tag, " saved method failures (", n_failures, " recorded entries).")
  )
  invisible(TRUE)
}

## Checks of the archived results.
lscar_unpack_methods <- function(z) {
  if (!is.null(z$method_results)) unname(z$method_results) else list(z)
}

lscar_validate_method <- function(v, r, group, expected_y = NULL) {
  m <- v$metrics
  d <- v$node_scores
  tag <- paste(group, r, if (is.data.frame(m)) m$method[1L] else "unknown")
  check <- function(ok, text) lscar_assert(ok, paste(tag, text))
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
    "invalid SCAR counts/labels"
  )
  check(
    identical(observed, c(rep(0L, 370), rep(1L, 30))) &&
      identical(hidden, 356:370),
    "linear SCAR uses the fixed hidden set 356:370"
  )
  if (!is.null(expected_y)) {
    check(
      identical(observed, as.integer(expected_y)),
      "observed labels differ from the common linear SCAR dataset"
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
    Overall_AUC = lscar_auc(score, truth),
    Hidden_AUC = lscar_auc(score[unlabelled], truth[unlabelled]),
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
      lscar_near(point, unlist(m[1L, names(point)]), 1e-8),
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
      lscar_near(calculated, unlist(m[1L, cols]), 1e-8),
      "coverage/interval scores disagree with saved node limits"
    )
  }
  m$Method <- lscar_canonical(m$method)
  m$group <- group
  d$Method <- lscar_canonical(d$method)
  d$group <- group
  list(metrics = m, scores = d)
}

# The labels of the linear SCAR design are fixed.
lscar_validate_dgp <- function(signature, tag) {
  cov <- signature$covariate_dgp
  if (is.null(cov) && !is.null(signature$dgp)) {
    cov <- signature$dgp$covariates
  }
  lscar_assert(
    is.list(cov) &&
      lscar_near(
        unlist(cov[c(
          "n_covariates",
          "background_mean",
          "signal_mean",
          "common_sd"
        )]),
        c(10, 0, 0.5, 1),
        0
      ),
    paste(tag, "wrong Gaussian data-generating specification")
  )
  net <- signature$network_dgp
  if (is.null(net) && !is.null(signature$dgp)) {
    net <- signature$dgp$network
  }
  if (!is.null(net)) {
    lscar_assert(
      lscar_near(
        unlist(net[c("p_00", "p_01", "p_11", "triad_closure_prob")]),
        c(0.10, 0.08, 0.25, 0),
        0
      ),
      paste(tag, "wrong SBM specification")
    )
  }
  invisible(TRUE)
}

lscar_read_group <- function(archive, group, expected_labels, reps = 1:100) {
  dir <- file.path(archive, group)
  sp <- lscar_group_spec()
  sp <- sp[sp$group == group, , drop = FALSE]
  paths <- file.path(dir, sp$checkpoint_dir, sprintf(sp$pattern, reps))
  expected <- lscar_registry()
  expected <- expected$Method[expected$group == group]
  out <- vector("list", length(reps))
  sig <- NULL
  for (j in seq_along(reps)) {
    r <- reps[j]
    x <- lscar_read(paths[j])
    lscar_assert(
      isTRUE(x$complete) && identical(as.integer(x$replicate), as.integer(r)),
      paste(group, r, "incomplete checkpoint")
    )
    if (is.null(sig)) {
      sig <- x$settings_signature
    }
    lscar_assert(
      identical(x$settings_signature, sig),
      paste(group, r, "different settings")
    )
    lscar_validate_dgp(sig, paste(group, r))
    lscar_assert(
      sig$master_seed == 20260718L && sig$n_rep == 100L,
      paste(group, r, "not the original experiment")
    )
    geo <- x$geometry_summary
    expected_seeds <- c(
      covariate_seed = 20260718 + 100000 + 1000 * r + 11,
      network_seed = 20260718 + 200000 + 1000 * r + 29,
      modularity_seed = 20260718 + 300000 + 1000 * r + 43
    )
    present <- intersect(names(expected_seeds), names(geo))
    lscar_assert(
      "covariate_seed" %in%
        present &&
        lscar_near(
          unlist(geo[1L, present, drop = FALSE]),
          expected_seeds[present],
          0
        ),
      paste(group, r, "different seed plan")
    )
    if (!group %in% c("Competitors", "nnPU")) {
      lscar_assert(
        lscar_near(
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
        lscar_unpack_methods(x),
        function(v) as.numeric(v$n_saved),
        numeric(1)
      )
      lscar_assert(all(saved == 2000L), paste(group, r, "retained draw count"))
    }
    if (group == "Bayesian_linear") {
      ev <- sig$evaluation$predictive_label_draws
      lscar_assert(
        is.character(ev) &&
          length(ev) == 1L &&
          grepl("Bernoulli(plogis", ev, fixed = TRUE),
        "Unexpected evaluation version of the Bayesian linear models."
      )
      ap <- sig$priors
      lp <- lscar_linear_priors()
      lscar_validate_linear_priors(ap, "Bayesian linear archive")
    }
    pp <- lapply(
      lscar_unpack_methods(x),
      lscar_validate_method,
      r = r,
      group = group,
      expected_y = expected_labels[[r]]
    )
    got <- unlist(lapply(pp, function(v) v$metrics$Method), use.names = FALSE)
    lscar_assert(
      !anyDuplicated(got) && setequal(got, expected),
      paste(group, r, "missing/extra methods")
    )
    out[[j]] <- pp
    lscar_assert_empty_method_failures(
      x[["method_failures", exact = TRUE]],
      paste(group, r)
    )
    if (!is.null(x$failure)) {
      stop(paste(group, r, "saved failure"), call. = FALSE)
    }
  }
  flat <- do.call(c, out)
  m <- lscar_bind(lapply(flat, `[[`, "metrics"))
  d <- lscar_bind(lapply(flat, `[[`, "scores"))
  csv <- utils::read.csv(
    file.path(dir, paste0(sp$prefix, "_metrics_by_replicate.csv")),
    stringsAsFactors = FALSE
  )
  for (method in unique(m$method)) {
    aa <- m[m$method == method, , drop = FALSE]
    bb <- csv[csv$method == method, , drop = FALSE]
    aa <- aa[order(aa$replicate), , drop = FALSE]
    bb <- bb[order(bb$replicate), , drop = FALSE]
    lscar_assert(
      nrow(aa) == 100L &&
        nrow(bb) == 100L &&
        !anyDuplicated(bb$replicate) &&
        setequal(bb$replicate, 1:100),
      paste(group, method, "CSV count")
    )
    for (metric in lscar_metrics()) {
      lscar_assert(
        lscar_near(aa[[metric]], bb[[metric]], 1e-8),
        paste(group, method, metric, "CSV/checkpoint mismatch")
      )
    }
  }
  cat(
    "PASS: ",
    group,
    ": 100 replications, ",
    length(expected),
    if (length(expected) == 1L) " method" else " methods",
    ", fixed linear-SCAR labels.\n",
    sep = ""
  )
  list(metrics = m, scores = d, signature = sig)
}

lscar_archive_manifest_check <- function(root, archive) {
  f <- utils::read.csv(
    file.path(lscar_source_root(root), "reference/ARCHIVE_MANIFEST.csv"),
    stringsAsFactors = FALSE
  )
  paths <- file.path(archive, f$file)
  missing <- !file.exists(paths)
  lscar_assert(
    !any(missing),
    paste(
      "Missing archive files:",
      paste(head(f$file[missing], 8L), collapse = "\n")
    )
  )
  lscar_assert(
    identical(as.character(lscar_md5(paths)), as.character(f$md5)),
    "Archive files differ from ARCHIVE_MANIFEST.csv."
  )
  cat("PASS: ", nrow(f), " archive files match the manifest.\n", sep = "")
  inventory <- utils::read.csv(
    file.path(
      root,
      "simulation_linear_SCAR_all18/reference/ARCHIVE_MANIFEST.csv"
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
  lscar_assert(
    setequal(actual_linear, expected_linear) &&
      length(actual_linear) == length(expected_linear),
    "Unexpected extra or missing Bayesian archive files"
  )
  invisible(TRUE)
}

lscar_collect_archive <- function(root, archive) {
  lscar_linear_archive_checks(archive)

  lscar_archive_manifest_check(root, archive)
  lscar_pilot_archive_checks(archive)
  labels <- vector("list", 100L)
  for (r in 1:100) {
    fn <- sprintf("replicate_%03d_geometry.rds", r)
    a <- lscar_read(file.path(archive, "Competitors/generated_geometries", fn))
    b <- lscar_read(file.path(archive, "nnPU/generated_geometries", fn))
    for (field in c(
      "X_raw",
      "X_std",
      "A_network",
      "Z_network_modularity",
      "T",
      "Y"
    )) {
      lscar_assert(
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
  groups <- lapply(lscar_groups(), function(g) {
    lscar_read_group(archive, g, labels)
  })
  names(groups) <- lscar_groups()
  metrics <- lscar_bind(lapply(groups, `[[`, "metrics"))
  scores <- lscar_bind(lapply(groups, `[[`, "scores"))
  lscar_assert(
    nrow(metrics) == 1800L && nrow(scores) == 720000L,
    "Incomplete all-method archive."
  )
  list(
    metrics = metrics,
    scores = scores,
    mode = "archived_full",
    replications = 100L
  )
}

lscar_pilot_archive_checks <- function(archive) {
  specs <- lscar_group_spec()
  total <- 0L
  for (group in setdiff(lscar_groups(), c("Competitors", "nnPU"))) {
    folder <- file.path(archive, group)
    files <- list.files(
      file.path(folder, "pilot_checkpoints"),
      pattern = "[.]rds$",
      full.names = TRUE
    )
    expected <- if (group == "Bayesian_linear") 20L else 5L
    lscar_assert(
      length(files) == expected,
      paste(group, "pilot checkpoint count")
    )
    for (path in files) {
      x <- lscar_read(path)
      sig <- x$settings_signature
      lscar_assert(
        is.list(x$tail) &&
          is.list(x$final_state) &&
          (group == "Bayesian_linear" || isTRUE(x$complete)),
        paste("Incomplete pilot", path)
      )
      lscar_assert(
        sig$pilot_burnin == 40000L && sig$pilot_tail_length == 5000L,
        paste("Pilot settings mismatch", path)
      )
      lscar_validate_dgp(sig, paste("pilot", group, basename(path)))
      if (group == "Bayesian_linear") {
        ap <- sig$priors
        lscar_validate_linear_priors(ap, "Bayesian linear archive")
      }
    }
    sp <- specs[specs$group == group, , drop = FALSE]
    stem <- sp$calibration_stem
    hash <- lscar_md5(file.path(folder, paste0(stem, ".rds")))
    first <- lscar_read(file.path(
      folder,
      sp$checkpoint_dir,
      sprintf(sp$pattern, 1L)
    ))
    lscar_assert(
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
