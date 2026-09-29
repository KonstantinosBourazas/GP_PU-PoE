## Tables of the prior sensitivity analysis over 100 datasets, rebuilt from the
## summaries of the 1,300 archived fits.

mc100_ids <- function() {
  c(
    "DEFAULT",
    "SIGMA_CONSERVATIVE",
    "SIGMA_DIFFUSE",
    "ELL_CONSERVATIVE",
    "ELL_DIFFUSE",
    "GAMMA_CONSERVATIVE",
    "GAMMA_DIFFUSE",
    "BETA0_CONSERVATIVE",
    "BETA0_DIFFUSE",
    "ETA_CONSERVATIVE",
    "ETA_DIFFUSE",
    "JOINT_CONSERVATIVE",
    "JOINT_DIFFUSE"
  )
}

mc100_groups <- function() {
  c("true zero", "hidden positive", "observed positive")
}

mc100_bind <- function(x) {
  z <- do.call(rbind, x)
  rownames(z) <- NULL
  z
}

mc100_assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

mc100_near <- function(x, y, tolerance = 1e-10) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  length(x) == length(y) &&
    all(is.finite(x)) &&
    all(is.finite(y)) &&
    all(abs(x - y) <= tolerance)
}

mc100_write <- function(object, path, csv = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".writing_", tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  if (csv) {
    utils::write.csv(object, tmp, row.names = FALSE)
  } else {
    saveRDS(object, tmp)
  }
  if (file.exists(path) && !file.remove(path)) {
    stop("Cannot replace output: ", path)
  }
  if (!file.rename(tmp, path)) {
    stop("Cannot finish writing: ", path)
  }
  invisible(path)
}

mc100_table <- function(x, path) mc100_write(x, path, csv = TRUE)

mc100_path <- function(archive, id, r) {
  file.path(archive, "checkpoints", id, sprintf("rep_%03d.rds", as.integer(r)))
}

mc100_read <- function(path) {
  tryCatch(readRDS(path), error = function(e) {
    stop(
      "Cannot read archived checkpoint: ",
      path,
      "\n",
      conditionMessage(e),
      call. = FALSE
    )
  })
}

mc100_manifest_check <- function(root, archive) {
  m <- utils::read.csv(
    file.path(root, "sensitivity_mc100", "reference", "ARCHIVE_MANIFEST.csv"),
    stringsAsFactors = FALSE
  )
  paths <- file.path(archive, m$file)
  missing <- !file.exists(paths)
  if (any(missing)) {
    stop(
      "Missing archive file(s):\n",
      paste(head(m$file[missing], 15L), collapse = "\n"),
      "\nThe archive folder must be complete and unchanged.",
      call. = FALSE
    )
  }
  md5 <- unname(tools::md5sum(paths))
  bad <- is.na(md5) | md5 != m$md5 | file.info(paths)$size != m$bytes
  if (any(bad)) {
    stop(
      "Archive files differ from ARCHIVE_MANIFEST.csv:\n",
      paste(head(m$file[bad], 15L), collapse = "\n"),
      call. = FALSE
    )
  }
  actual <- list.files(
    file.path(archive, "checkpoints"),
    pattern = "[.]rds$",
    recursive = TRUE,
    full.names = FALSE
  )
  expected <- sub(
    "^checkpoints/",
    "",
    m$file[startsWith(m$file, "checkpoints/")]
  )
  actual <- gsub("\\\\", "/", actual)
  mc100_assert(
    length(actual) == 1300L && setequal(actual, expected),
    "Expected the 1,300 checkpoint files of the manifest; found a different set."
  )
  m
}

mc100_auc <- function(p, truth) {
  np <- sum(truth == 1L)
  nn <- sum(truth == 0L)
  mc100_assert(np > 0L && nn > 0L, "AUC requires both classes.")
  (sum(rank(p, ties.method = "average")[truth == 1L]) - np * (np + 1) / 2) /
    (np * nn)
}

mc100_eta_mean_row <- function(d, id) {
  data.frame(
    Setting = id,
    Replications = nrow(d),
    True_eta = 1 / 3,
    Mean_posterior_mean = mean(d$eta_posterior_mean),
    Mean_HPD80_lower = mean(d$eta_HPD80_lower),
    Mean_HPD80_upper = mean(d$eta_HPD80_upper),
    Mean_HPD80_width = mean(d$eta_HPD80_width),
    Coverage_HPD80 = mean(d$HPD80_covers_true),
    stringsAsFactors = FALSE
  )
}

mc100_extra_eta <- function(d, id) {
  estimate <- d$eta_posterior_mean
  truth <- 1 / 3
  coverage <- mean(d$HPD80_covers_true)
  data.frame(
    Setting = id,
    Replications = nrow(d),
    Bias_posterior_mean = mean(estimate - truth),
    RMSE_posterior_mean = sqrt(mean((estimate - truth)^2)),
    SD_across_posterior_means = stats::sd(estimate),
    MCSE_average_posterior_mean = stats::sd(estimate) / sqrt(nrow(d)),
    Coverage_HPD80 = coverage,
    MCSE_coverage = sqrt(coverage * (1 - coverage) / nrow(d)),
    stringsAsFactors = FALSE
  )
}

mc100_check_fit <- function(z, id, r, experiment, default_hash) {
  tag <- paste0(id, "/", sprintf("rep_%03d", r))
  ck <- function(ok, message) mc100_assert(ok, paste0(tag, ": ", message))
  required <- c(
    "complete",
    "experiment",
    "fit_id",
    "replicate",
    "toy_hash",
    "fit_signature",
    "fit_spec",
    "eta_prior",
    "beta0_prior",
    "node_summary",
    "eta_summary",
    "diagnostics",
    "prior_calibration",
    "parameter_summary",
    "core_parameter_ess",
    "node_probability_ess",
    "node_ess_summary"
  )
  ck(is.list(z) && all(required %in% names(z)), "Incomplete checkpoint schema.")
  ck(
    isTRUE(z$complete) &&
      identical(z$fit_id, id) &&
      identical(z$replicate, as.integer(r)),
    "Identity/completion mismatch."
  )
  ck(
    identical(z$experiment, experiment),
    "Different archived experiment metadata."
  )
  ck(
    identical(z$toy_hash, default_hash),
    "Dataset hash differs from DEFAULT in this replication."
  )
  s <- z$fit_signature
  sp <- experiment$seed_plan[r, ]
  ck(identical(s$toy_data_hash, z$toy_hash), "Internal dataset hash mismatch.")
  ck(
    s$toy_seed == sp$data_seed &&
      s$modularity_seed == sp$modularity_seed &&
      z$fit_spec$mcmc_seed == sp$mcmc_seed,
    "Seed-plan mismatch."
  )
  expected <- experiment$specifications[[id]]
  expected$mcmc_seed <- as.integer(sp$mcmc_seed)
  ck(
    identical(z$fit_spec, expected) && identical(s$fit_spec, z$fit_spec),
    "Prior specification mismatch."
  )
  ck(
    identical(s$priors$eta, z$eta_prior) &&
      identical(s$priors$beta0, z$beta0_prior),
    "Prior metadata mismatch."
  )
  ck(
    mc100_near(
      unlist(s$MCMC[c(
        "blocked_burnin",
        "joint_burnin",
        "posterior_iterations",
        "thin",
        "retained"
      )]),
      c(10000, 20000, 150000, 30, 5000),
      0
    ),
    "Not the archived full schedule."
  )
  ck(
    identical(s$requested_interval$type, "HPD") &&
      s$requested_interval$probability == 0.80,
    "Not the 80% HPD summaries."
  )
  d <- z$node_summary
  dg <- z$diagnostics
  et <- z$eta_summary
  ck(
    is.data.frame(d) &&
      nrow(d) == 100L &&
      !anyDuplicated(d$node) &&
      setequal(d$node, 1:100),
    "Bad node index."
  )
  d <- d[match(1:100, d$node), , drop = FALSE]
  ck(
    all(d$T == c(rep(0L, 70), rep(1L, 30))) &&
      all(d$Y == c(rep(0L, 80), rep(1L, 20))) &&
      identical(as.character(d$node_group), rep(mc100_groups(), c(70, 10, 20))),
    "Different class design."
  )
  for (prefix in c("latent", "observed", "conditional_T")) {
    p <- d[[paste0(prefix, "_probability_mean")]]
    lo <- d[[paste0(prefix, "_HPD80_lower")]]
    hi <- d[[paste0(prefix, "_HPD80_upper")]]
    w <- d[[paste0(prefix, "_HPD80_width")]]
    ck(
      length(p) == 100L &&
        length(lo) == 100L &&
        length(hi) == 100L &&
        all(is.finite(c(p, lo, hi, w))) &&
        all(p >= 0 & p <= 1) &&
        all(lo >= 0 & hi <= 1 & lo <= hi) &&
        mc100_near(w, hi - lo),
      "Invalid probability/HPD summaries."
    )
  }
  ck(
    mc100_near(d$posterior_predictive_mean, d$latent_probability_mean, 0) &&
      mc100_near(d$HPD80_lower, d$latent_HPD80_lower, 0) &&
      mc100_near(d$HPD80_upper, d$latent_HPD80_upper, 0),
    "Figure target alias mismatch."
  )
  ck(
    nrow(et) == 1L &&
      et$eta_true == 1 / 3 &&
      all(is.finite(unlist(et[c(
        "eta_posterior_mean",
        "eta_posterior_sd",
        "eta_HPD80_lower",
        "eta_HPD80_upper",
        "eta_HPD80_width"
      )]))),
    "Invalid eta row."
  )
  ck(
    et$eta_HPD80_lower >= 0 &&
      et$eta_HPD80_upper <= 1 &&
      et$eta_HPD80_lower <= et$eta_HPD80_upper &&
      mc100_near(et$eta_HPD80_width, et$eta_HPD80_upper - et$eta_HPD80_lower) &&
      et$HPD80_covers_true ==
        as.integer(et$eta_HPD80_lower <= 1 / 3 && et$eta_HPD80_upper >= 1 / 3),
    "Eta HPD width/coverage mismatch."
  )
  ep <- z$parameter_summary[
    z$parameter_summary$parameter == "eta",
    ,
    drop = FALSE
  ]
  ck(
    nrow(ep) == 1L &&
      mc100_near(
        unlist(et[c(
          "eta_posterior_mean",
          "eta_posterior_sd",
          "eta_HPD80_lower",
          "eta_HPD80_upper"
        )]),
        unlist(ep[c(
          "posterior_mean",
          "posterior_sd",
          "HPD80_lower",
          "HPD80_upper"
        )])
      ),
    "Eta cross-summary mismatch."
  )
  ck(
    nrow(dg) == 1L &&
      mc100_near(
        unlist(dg[c(
          "blocked_burnin",
          "joint_burnin",
          "posterior_iterations",
          "thinning",
          "retained_draws"
        )]),
        c(10000, 20000, 150000, 30, 5000),
        0
      ),
    "Diagnostics schedule mismatch."
  )
  p <- d$latent_probability_mean
  clipped <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  metrics <- c(
    mc100_auc(p, d$T),
    mc100_auc(p[1:80], d$T[1:80]),
    -mean(d$T * log(clipped) + (1 - d$T) * log(1 - clipped)),
    mean((p - d$T)^2),
    mean(d$latent_HPD80_width)
  )
  ck(
    mc100_near(
      metrics,
      unlist(dg[c(
        "overall_AUC_latent_mean",
        "hidden_AUC_latent_mean",
        "overall_LogLoss_latent_mean",
        "overall_Brier_latent_mean",
        "mean_HPD80_width"
      )])
    ),
    "Metrics do not match archived node means."
  )
  ce <- z$core_parameter_ess
  ne <- z$node_probability_ess
  ck(
    nrow(ce) == 8L &&
      nrow(ne) == 100L &&
      all(ce$retained_draws == 5000L) &&
      all(ne$retained_draws == 5000L) &&
      all(is.finite(c(ce$ESS, ne$ESS))) &&
      all(c(ce$ESS, ne$ESS) > 0),
    "Invalid archived ESS/counts."
  )
  ne_summary <- c(
    min(ne$ESS),
    unname(stats::quantile(ne$ESS, .25)),
    stats::median(ne$ESS),
    mean(ne$ESS),
    unname(stats::quantile(ne$ESS, .75)),
    max(ne$ESS)
  )
  ck(
    mc100_near(
      ne_summary,
      unlist(z$node_ess_summary[c(
        "min_ESS",
        "q25_ESS",
        "median_ESS",
        "mean_ESS",
        "q75_ESS",
        "max_ESS"
      )]),
      1e-8
    ),
    "Nodewise ESS summary mismatch."
  )
  invisible(TRUE)
}

mc100_compare_table <- function(actual, expected, label, tolerance = 1e-10) {
  mc100_assert(
    nrow(actual) == nrow(expected) && identical(names(actual), names(expected)),
    paste0(label, ": schema mismatch.")
  )
  max_error <- 0
  for (key in names(expected)) {
    if (is.numeric(expected[[key]])) {
      mc100_assert(
        mc100_near(actual[[key]], expected[[key]], tolerance),
        paste0(label, ": mismatch in ", key)
      )
      max_error <- max(
        max_error,
        abs(as.numeric(actual[[key]]) - expected[[key]])
      )
    } else {
      mc100_assert(
        identical(as.character(actual[[key]]), as.character(expected[[key]])),
        paste0(label, ": mismatch in ", key)
      )
    }
  }
  data.frame(
    check = label,
    passed = TRUE,
    max_absolute_error = max_error,
    stringsAsFactors = FALSE
  )
}

mc100_combine <- function(root, archive, out) {
  ids <- mc100_ids()
  reps <- 1:100
  cat("Checking the 1,300 archived checkpoints and 39 archived tables...\n")
  manifest <- mc100_manifest_check(root, archive)
  first <- mc100_read(mc100_path(archive, "DEFAULT", 1L))
  experiment <- first$experiment
  mc100_assert(
    identical(experiment$driver_version, "MC100_HPD80_v1"),
    "Unexpected archived driver version."
  )
  mc100_assert(
    identical(names(experiment$specifications), ids),
    "Unexpected archived prior settings/order."
  )
  reference <- file.path(
    root,
    "sensitivity_mc100",
    "reference",
    "prior_sensitivity_single_run_original.R"
  )
  mc100_assert(
    identical(unname(tools::md5sum(reference)), experiment$source_md5),
    "The MD5 of reference/prior_sensitivity_single_run_original.R differs from the one recorded in the archive."
  )
  seed_plan <- experiment$seed_plan
  mc100_assert(
    nrow(seed_plan) == 100L &&
      mc100_near(seed_plan$replicate, reps, 0) &&
      mc100_near(seed_plan$data_seed, 53L + 100003L * (0:99), 0) &&
      mc100_near(seed_plan$modularity_seed, 53043L + 100003L * (0:99), 0) &&
      mc100_near(seed_plan$mcmc_seed, 530401L + 100003L * (0:99), 0),
    "Archived seed plan mismatch."
  )
  defaults <- lapply(reps, function(r) {
    mc100_read(mc100_path(archive, "DEFAULT", r))
  })
  default_hashes <- vapply(defaults, function(z) z$toy_hash, character(1))
  mc100_assert(
    length(unique(default_hashes)) == 100L,
    "Expected 100 distinct archived dataset hashes."
  )
  ## Units keep their indices in the simulation, without sorting by posterior probability.
  orders <- lapply(reps, function(r) 1:100)
  node_positions <- mc100_bind(lapply(reps, function(r) {
    data.frame(
      replicate = r,
      plot_position = 1:100,
      original_node = orders[[r]],
      node_group = rep(mc100_groups(), c(70, 10, 20)),
      stringsAsFactors = FALSE
    )
  }))
  prediction <- eta_rows <- eta_rep <- eta_extra <- robustness <- topk <- inventory <- checks <- list()
  core_ess <- node_ess <- node_variation <- list()
  table_names <- c(
    "node_summary",
    "eta_summary",
    "diagnostics",
    "prior_calibration",
    "parameter_summary",
    "core_parameter_ess",
    "node_probability_ess",
    "node_ess_summary"
  )
  for (id in ids) {
    cat("Reading ", id, " (100 full-run summaries)...\n", sep = "")
    zz <- if (id == "DEFAULT") {
      defaults
    } else {
      lapply(reps, function(r) mc100_read(mc100_path(archive, id, r)))
    }
    for (r in reps) {
      mc100_check_fit(zz[[r]], id, r, experiment, default_hashes[r])
    }
    folder <- file.path(out, "tables", id)
    for (key in table_names) {
      mc100_table(
        mc100_bind(lapply(zz, `[[`, key)),
        file.path(folder, paste0(key, "_by_replication.csv"))
      )
    }
    et <- mc100_bind(lapply(zz, `[[`, "eta_summary"))
    eta_rep[[id]] <- et
    eta_rows[[id]] <- mc100_eta_mean_row(et, id)
    eta_extra[[id]] <- mc100_extra_eta(et, id)
    old_per_rep <- utils::read.csv(
      file.path(archive, "tables", id, "eta_by_replication.csv"),
      stringsAsFactors = FALSE
    )
    old_per_rep <- old_per_rep[order(old_per_rep$replicate), , drop = FALSE]
    checks[[length(checks) + 1L]] <- mc100_compare_table(
      et,
      old_per_rep,
      paste0(id, ": archived eta table by replication")
    )
    old_mean <- utils::read.csv(
      file.path(archive, "tables", id, "eta_means_full_precision.csv"),
      stringsAsFactors = FALSE
    )
    checks[[length(checks) + 1L]] <- mc100_compare_table(
      eta_rows[[id]],
      old_mean,
      paste0(id, ": archived table of eta means")
    )
    display <- eta_rows[[id]][c(
      "Setting",
      "Replications",
      "Mean_posterior_mean",
      "Mean_HPD80_lower",
      "Mean_HPD80_upper"
    )]
    for (key in names(display)[-(1:2)]) {
      display[[key]] <- sprintf("%.2f", display[[key]])
    }
    old_display <- utils::read.csv(
      file.path(archive, "tables", id, "eta_means_2dp.csv"),
      stringsAsFactors = FALSE,
      colClasses = "character"
    )
    checks[[length(checks) + 1L]] <- mc100_compare_table(
      display,
      old_display,
      paste0(id, ": archived table of rounded eta means")
    )
    mc100_table(
      eta_rows[[id]],
      file.path(folder, "eta_means_full_precision.csv")
    )
    mc100_table(display, file.path(folder, "eta_means_2dp.csv"))
    aligned <- lapply(reps, function(r) {
      d <- zz[[r]]$node_summary
      as.matrix(d[
        match(orders[[r]], d$node),
        c(
          "latent_probability_mean",
          "latent_HPD80_lower",
          "latent_HPD80_upper"
        ),
        drop = FALSE
      ])
    })
    av <- Reduce(`+`, aligned) / 100
    prediction[[id]] <- data.frame(
      fit_id = id,
      position = 1:100,
      node_group = rep(mc100_groups(), c(70, 10, 20)),
      n_replications = 100L,
      Mean_posterior_mean = av[, 1L],
      Mean_HPD80_lower = av[, 2L],
      Mean_HPD80_upper = av[, 3L],
      stringsAsFactors = FALSE
    )
    mc100_table(
      prediction[[id]],
      file.path(folder, "predictive_means_original_index.csv")
    )
    literal <- Reduce(
      `+`,
      lapply(zz, function(z) {
        as.matrix(z$node_summary[
          match(1:100, z$node_summary$node),
          c(
            "latent_probability_mean",
            "latent_HPD80_lower",
            "latent_HPD80_upper"
          ),
          drop = FALSE
        ])
      })
    ) /
      100
    literal_df <- prediction[[id]]
    literal_df[, c(
      "Mean_posterior_mean",
      "Mean_HPD80_lower",
      "Mean_HPD80_upper"
    )] <- literal
    names(literal_df)[names(literal_df) == "position"] <- "original_node"
    mc100_table(
      literal_df,
      file.path(folder, "predictive_means_by_original_node.csv")
    )
    rob <- tp <- list()
    for (r in reps) {
      d <- zz[[r]]$node_summary
      b <- defaults[[r]]$node_summary
      p <- d$latent_probability_mean[match(1:100, d$node)]
      p0 <- b$latent_probability_mean[match(1:100, b$node)]
      dg <- zz[[r]]$diagnostics
      rob[[r]] <- mc100_robustness_row(zz[[r]], defaults[[r]], id, r)
      tp[[r]] <- mc100_bind(lapply(c(10L, 20L, 30L), function(k) {
        overlap <- length(intersect(
          order(-p, 1:100)[seq_len(k)],
          order(-p0, 1:100)[seq_len(k)]
        ))
        data.frame(
          fit_id = id,
          replicate = r,
          k = k,
          overlap_count = overlap,
          overlap_fraction = overlap / k
        )
      }))
    }
    robustness[[id]] <- mc100_bind(rob)
    topk[[id]] <- mc100_bind(tp)
    mc100_table(
      robustness[[id]],
      file.path(folder, "robustness_by_replication.csv")
    )
    mc100_table(
      topk[[id]],
      file.path(folder, "topk_overlap_by_replication.csv")
    )
    keys <- setdiff(names(robustness[[id]]), c("fit_id", "replicate"))
    mc100_table(
      data.frame(
        Metric = keys,
        Mean = vapply(robustness[[id]][keys], mean, numeric(1)),
        SD = vapply(robustness[[id]][keys], stats::sd, numeric(1))
      ),
      file.path(folder, "robustness_MC_mean_sd.csv")
    )
    node_variation[[id]] <- mc100_original_node_variation(zz, id)
    core_ess[[id]] <- mc100_bind(lapply(zz, `[[`, "core_parameter_ess"))
    node_ess[[id]] <- mc100_bind(lapply(zz, `[[`, "node_ess_summary"))
    inventory[[id]] <- data.frame(
      fit_id = id,
      n_replications = length(zz),
      retained_draws_per_fit = 5000L,
      first_completed = min(vapply(
        zz,
        function(z) z$completed_at,
        character(1)
      )),
      last_completed = max(vapply(
        zz,
        function(z) z$completed_at,
        character(1)
      )),
      stringsAsFactors = FALSE
    )
  }
  reports <- mc100_auc_reports(mc100_bind(robustness))
  combined <- list(
    format_version = "mc100_original_node_order",
    experiment = experiment,
    prediction_index = "original_node_no_probability_ordering",
    provenance = list(
      origin = "Summaries of the 1,300 archived full fits",
      source_reference_md5 = unname(tools::md5sum(reference)),
      plotting_session_info = capture.output(sessionInfo()),
      created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ),
    input_manifest = manifest,
    plot_means = mc100_bind(prediction),
    eta_means = mc100_bind(eta_rows),
    eta_by_replication = mc100_bind(eta_rep),
    eta_extra_metrics = mc100_bind(eta_extra),
    node_positions = node_positions,
    robustness = mc100_bind(robustness),
    topk_overlap = mc100_bind(topk),
    auc_reports = reports,
    node_variation = mc100_bind(node_variation),
    core_parameter_ess = mc100_bind(core_ess),
    node_ess_summary = mc100_bind(node_ess),
    interval_note = paste(
      "HPD limits were calculated within each fit and then averaged across 100 datasets.",
      "These averaged limits are not a pooled posterior interval or a confidence interval."
    )
  )
  for (item in list(
    c("plot_means", "expected_prediction_means.csv"),
    c("eta_means", "expected_eta_means.csv")
  )) {
    expected <- utils::read.csv(
      file.path(root, "sensitivity_mc100", "reference", item[2L]),
      stringsAsFactors = FALSE
    )
    checks[[length(checks) + 1L]] <- mc100_compare_table(
      combined[[item[1L]]],
      expected,
      paste0("Reference table: ", item[1L])
    )
  }
  for (item in list(
    c("summary", "expected_auc_summary.csv"),
    c("pairwise_correlations", "expected_auc_pairwise.csv"),
    c("range_summary", "expected_auc_range_summary.csv"),
    c("ranking_stability", "expected_ranking_stability.csv")
  )) {
    expected <- utils::read.csv(
      file.path(root, "sensitivity_mc100", "reference", item[2L]),
      stringsAsFactors = FALSE
    )
    checks[[length(checks) + 1L]] <- mc100_compare_table(
      reports[[item[1L]]],
      expected,
      paste0("Reference table: ", item[1L]),
      tolerance = 1e-9
    )
  }
  mc100_table(
    combined$node_variation,
    file.path(
      out,
      "combined",
      "MC100_posterior_mean_variation_by_original_node.csv"
    )
  )
  mc100_save_auc_reports(reports, file.path(out, "combined"))
  mc100_table(
    mc100_bind(checks),
    file.path(out, "verification", "comparison_results.csv")
  )
  mc100_table(
    mc100_bind(inventory),
    file.path(out, "verification", "archive_completeness.csv")
  )
  mc100_table(manifest, file.path(out, "verification", "input_manifest.csv"))
  mc100_table(seed_plan, file.path(out, "combined", "seed_plan.csv"))
  mc100_table(
    node_positions,
    file.path(out, "combined", "original_node_positions_by_replication.csv")
  )
  mc100_table(
    combined$plot_means,
    file.path(out, "combined", "ALL_predictive_mean_HPD80.csv")
  )
  mc100_table(
    combined$eta_means,
    file.path(out, "combined", "ALL_eta_means_full_precision.csv")
  )
  mc100_table(
    combined$eta_extra_metrics,
    file.path(out, "combined", "ALL_eta_bias_RMSE_MCSE.csv")
  )
  mc100_table(
    combined$eta_by_replication,
    file.path(out, "combined", "ALL_eta_by_replication.csv")
  )
  mc100_table(
    combined$robustness,
    file.path(out, "combined", "ALL_robustness_by_replication.csv")
  )
  mc100_write(
    combined,
    file.path(out, "combined", "MC100_combined_summaries.rds")
  )
  cat(
    "Checked: 100 replications x 13 settings, metadata, common datasets and summaries.\n"
  )
  cat("The HPD intervals and ESS of each fit are the archived values.\n")
  print(combined$eta_means, row.names = FALSE)
  invisible(combined)
}
