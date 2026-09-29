## Checks of the archived results.
scar_unpack_methods <- function(z) {
  if (!is.null(z$method_results)) unname(z$method_results) else list(z)
}

scar_validate_method <- function(v, r, group) {
  m <- v$metrics
  d <- v$node_scores
  tag <- paste(group, r, if (is.data.frame(m)) m$method[1L] else "unknown")
  check <- function(x, s) scar_assert(x, paste(tag, s))
  check(
    is.data.frame(m) && nrow(m) == 1L && is.data.frame(d) && nrow(d) == 400L,
    "bad result dimensions"
  )
  check(
    all(m$replicate == r) &&
      all(d$replicate == r) &&
      !anyDuplicated(d$node) &&
      setequal(d$node, 1:400),
    "bad replication/node identities"
  )
  d <- d[match(1:400, d$node), , drop = FALSE]
  truth <- c(rep(0L, 355), rep(1L, 45))
  observed <- c(rep(0L, 370), rep(1L, 30))
  check(
    identical(as.integer(d$T), truth) && identical(as.integer(d$Y), observed),
    "different nonlinear SCAR truth"
  )
  p <- as.numeric(d$probability)
  score <- if ("decision_score" %in% names(d)) d$decision_score else p
  check(all(is.finite(c(p, score))) && all(p >= 0 & p <= 1), "invalid scores")
  rk <- rank(-score, ties.method = "average")
  ru <- rank(-score[1:370], ties.method = "average")
  clip <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  point <- c(
    Overall_AUC = scar_auc(score, truth),
    Hidden_AUC = scar_auc(score[1:370], truth[1:370]),
    Overall_Recall_at_45 = sum(truth == 1L & rk <= 45) / 45,
    Hidden_Recall_at_15 = sum(truth[1:370] == 1L & ru <= 15) / 15,
    Overall_LogLoss = -mean(truth * log(clip) + (1 - truth) * log(1 - clip)),
    Hidden_LogLoss = -mean(log(clip[356:370])),
    Overall_Brier = mean((p - truth)^2),
    Hidden_Brier = mean((p[356:370] - 1)^2)
  )
  check(
    all(names(point) %in% names(m)) &&
      scar_near(point, unlist(m[1, names(point)]), 1e-8),
    "point metrics do not match archived node scores"
  )
  for (level in c(80L, 95L)) {
    lo_name <- paste0("PI", level, "_lower")
    hi_name <- paste0("PI", level, "_upper")
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
      all(c(lo_name, hi_name) %in% names(d)) && all(cols %in% names(m)),
      "missing prediction-interval fields"
    )
    lo <- d[[lo_name]]
    hi <- d[[hi_name]]
    if (all(is.na(lo)) && all(is.na(hi))) {
      check(
        all(is.na(unlist(m[1, cols]))),
        "unavailable uncertainty must remain NA"
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
    coverage <- as.numeric(lo <= truth & truth <= hi)
    alpha <- 1 - level / 100
    iscore <- hi -
      lo +
      2 / alpha * (lo - truth) * (truth < lo) +
      2 / alpha * (truth - hi) * (truth > hi)
    calc <- c(
      mean(coverage),
      mean(coverage[356:370]),
      mean(iscore),
      mean(iscore[356:370])
    )
    check(
      scar_near(calc, unlist(m[1, cols]), 1e-8),
      "coverage/interval scores do not match stored node limits"
    )
  }
  m$Method <- scar_canonical(m$method)
  m$group <- group
  d$Method <- scar_canonical(d$method)
  d$group <- group
  list(metrics = m, scores = d)
}

scar_read_group <- function(archive, group, full = TRUE, reps = 1:100) {
  dir <- file.path(archive, group)
  spec <- scar_group_spec()
  sp <- spec[spec$group == group, , drop = FALSE]
  if (group == "Competitors" && full) {
    ## The complete results bundle holds the evaluation of the seven competitors.
    ## Their replicate_checkpoints are kept in the archive but not read.
    p <- file.path(
      dir,
      "SIMULATION_100_COMPETING_METHODS_COMPLETE_RESULTS.RData"
    )
    scar_assert(file.exists(p), paste("Missing competitor results bundle", p))
    e <- new.env(parent = emptyenv())
    load(p, envir = e)
    scar_assert(
      exists("replication_results", e, inherits = FALSE),
      "Competitor bundle lacks replication_results."
    )
    z <- e$replication_results
    scar_assert(
      length(z) == 100L,
      "Competitor bundle does not contain all 100 replications."
    )
    z <- z[reps]
  } else {
    paths <- file.path(dir, sp$checkpoint_dir, sprintf(sp$pattern, reps))
    z <- lapply(paths, scar_read)
  }
  expected <- scar_registry()
  expected <- expected$Method[expected$group == group]
  sig <- z[[1L]]$settings_signature
  out <- vector("list", length(reps))
  checks <- list()
  for (j in seq_along(reps)) {
    r <- reps[j]
    x <- z[[j]]
    scar_assert(
      isTRUE(x$complete) && identical(as.integer(x$replicate), as.integer(r)),
      paste(group, r, "incomplete checkpoint")
    )
    scar_assert(
      identical(x$settings_signature, sig),
      paste(group, r, "inconsistent settings signature")
    )
    geom <- x$geometry_summary
    expected_seeds <- c(
      covariate_seed = 20260718 + 100000 + 1000 * r + 11,
      network_seed = 20260718 + 200000 + 1000 * r + 29,
      modularity_seed = 20260718 + 300000 + 1000 * r + 43
    )
    present <- intersect(names(expected_seeds), names(geom))
    scar_assert(
      length(present) > 0L &&
        scar_near(
          unlist(geom[1, present, drop = FALSE]),
          expected_seeds[present],
          0
        ),
      paste(group, r, "seed metadata mismatch")
    )
    if (full && !group %in% c("Competitors", "nnPU")) {
      scar_assert(
        scar_near(
          unlist(sig[c(
            "production_burnin",
            "production_sampling",
            "production_thin"
          )]),
          c(15000, 40000, 20),
          0
        ),
        paste(group, "not archived full schedule")
      )
      retained <- vapply(
        scar_unpack_methods(x),
        function(v) as.numeric(v$n_saved),
        numeric(1)
      )
      scar_assert(
        all(retained == 2000),
        paste(group, r, "retained count mismatch")
      )
    }
    pp <- lapply(
      scar_unpack_methods(x),
      scar_validate_method,
      r = r,
      group = group
    )
    got <- unlist(lapply(pp, function(v) v$metrics$Method), use.names = FALSE)
    scar_assert(
      !anyDuplicated(got) && setequal(got, expected),
      paste(group, r, "missing/extra methods")
    )
    out[[j]] <- pp
    if (full && group == "Bayesian_linear") {
      ev <- x$settings_signature$evaluation
      scar_assert(
        !is.null(ev$predictive_label_draws) &&
          grepl("Bernoulli(plogis", ev$predictive_label_draws, fixed = TRUE),
        "Unexpected evaluation version of the Bayesian linear models."
      )
      lp <- scar_linear_priors()
      ap <- x$settings_signature$priors
      scar_validate_linear_priors(ap, "Bayesian linear archive")
    }
  }
  flat <- do.call(c, out)
  m <- scar_bind(lapply(flat, `[[`, "metrics"))
  d <- scar_bind(lapply(flat, `[[`, "scores"))
  if (full) {
    file <- file.path(dir, paste0(sp$prefix, "_metrics_by_replicate.csv"))
    csv <- utils::read.csv(file, stringsAsFactors = FALSE)
    for (method in unique(m$method)) {
      aa <- m[m$method == method, , drop = FALSE]
      bb <- csv[csv$method == method, , drop = FALSE]
      aa <- aa[order(aa$replicate), , drop = FALSE]
      bb <- bb[order(bb$replicate), , drop = FALSE]
      scar_assert(
        nrow(aa) == 100L && nrow(bb) == 100L && !anyDuplicated(bb$replicate),
        paste(group, method, "CSV count")
      )
      for (metric in scar_metrics()) {
        scar_assert(
          scar_near(aa[[metric]], bb[[metric]]),
          paste(group, method, metric, "CSV/checkpoint difference")
        )
      }
    }
  }
  cat(
    "PASS: ",
    group,
    ": ",
    length(reps),
    " replications, ",
    length(expected),
    if (length(expected) == 1L) " method.\n" else " methods.\n",
    sep = ""
  )
  list(
    metrics = m,
    scores = d,
    signature = sig,
    replications = length(reps),
    methods = length(expected)
  )
}

scar_archive_manifest_check <- function(root, archive) {
  f <- utils::read.csv(
    file.path(scar_source_root(root), "reference/ARCHIVE_MANIFEST.csv"),
    stringsAsFactors = FALSE
  )
  paths <- file.path(archive, f$file)
  missing <- !file.exists(paths)
  scar_assert(
    !any(missing),
    paste(
      "Missing archive files:",
      paste(head(f$file[missing], 8), collapse = "\n")
    )
  )
  actual <- scar_md5(paths)
  scar_assert(
    identical(as.character(actual), as.character(f$md5)),
    "Archive files differ from ARCHIVE_MANIFEST.csv."
  )
  cat("PASS: ", nrow(f), " archive files match the manifest.\n", sep = "")
  inventory <- utils::read.csv(
    file.path(
      root,
      "simulation_nonlinear_SCAR_all18/reference/ARCHIVE_MANIFEST.csv"
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
  scar_assert(
    setequal(actual_linear, expected_linear) &&
      length(actual_linear) == length(expected_linear),
    "Unexpected extra or missing Bayesian archive files"
  )
  invisible(TRUE)
}

scar_collect_archive <- function(root, archive) {
  scar_linear_archive_checks(archive)

  scar_archive_manifest_check(root, archive)
  for (r in 1:100) {
    a <- scar_read(file.path(
      archive,
      "Competitors/generated_geometries",
      sprintf("replicate_%03d_geometry.rds", r)
    ))
    b <- scar_read(file.path(
      archive,
      "nnPU/generated_geometries",
      sprintf("replicate_%03d_geometry.rds", r)
    ))
    for (field in c("X_std", "A_network", "T", "Y")) {
      scar_assert(
        isTRUE(all.equal(
          a[[field]],
          b[[field]],
          tolerance = 0,
          check.attributes = FALSE
        )),
        paste("Saved competitor/nnPU geometry differs", r, field)
      )
    }
  }
  cat(
    "PASS: 100 datasets shared by the competitors and nnPU; the GP archives keep seeds and summaries, not all X and A matrices.\n"
  )
  groups <- lapply(scar_groups(), function(g) scar_read_group(archive, g))
  names(groups) <- scar_groups()
  metrics <- scar_bind(lapply(groups, `[[`, "metrics"))
  scores <- scar_bind(lapply(groups, `[[`, "scores"))
  scar_assert(
    nrow(metrics) == 1800L && nrow(scores) == 720000L,
    "The all-method archive is incomplete."
  )
  list(
    metrics = metrics,
    scores = scores,
    groups = groups,
    mode = "archived_full",
    replications = 100L
  )
}
