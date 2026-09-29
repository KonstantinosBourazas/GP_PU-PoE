## Priors of the Bayesian linear models. The calibration itself is computed
## from reference/sources/Bayesian_linear.R.
lscar_linear_priors <- function() {
  list(
    beta0 = list(mean = -3.40, sd = .61),
    eta = list(a_eta = 2, b_eta = 6.39),
    tau_scale_rule = list(
      source = "internal geometry calibration, no application quantities",
      total_sd_center = 2,
      pair_variance_shares = c(cov = .5, network = .5),
      single_effective_sd_center = 2,
      halft_df = 3,
      mean_rms_95_mult = 4.16,
      round_active_priors = TRUE,
      prior_round_digits = 2L,
      b_cov_fixed = sqrt(10 * 399 / 400),
      tau_scale_cov_single = .83,
      network_scale = "recomputed per replication from the modularity dimension"
    ),
    nu_tau = 3
  )
}

lscar_validate_linear_priors <- function(priors, tag = "Bayesian linear") {
  expected <- lscar_linear_priors()
  lscar_assert(
    is.list(priors) && setequal(names(priors), names(expected)),
    paste(tag, "unrecognized prior schema")
  )
  for (n in names(expected)) {
    lscar_assert(
      isTRUE(all.equal(
        priors[[n]],
        expected[[n]],
        tolerance = 1e-12,
        check.attributes = TRUE
      )),
      paste(tag, "prior mismatch:", n)
    )
  }
  invisible(TRUE)
}

lscar_validate_linear_environment <- function(e) {
  lscar_validate_linear_priors(
    list(
      beta0 = e$beta0_prior_for_sampler,
      eta = e$eta_prior_for_sampler,
      tau_scale_rule = e$TAU_SCALE_RULE,
      nu_tau = e$NU_TAU
    ),
    "Imported source"
  )
  # Deterministic standardized matrices.
  i <- seq_len(400L)
  X <- scale(outer(i, seq_len(10L), function(i, j) {
    sin(i * j / 17) + cos(i / (j + 2))
  }))
  for (q in c(1L, 2L, 3L, 5L, 10L, 20L)) {
    Z <- scale(outer(i, seq_len(q), function(i, j) {
      cos(i * j / 23) + sin(i / (j + 3))
    }))
    specs <- e$make_model_specs(list(
      X_cov = X,
      Z_net = Z,
      X_cov_net = cbind(X, Z)
    ))
    lscar_assert(
      identical(names(specs), c("cov", "cov_pu", "cov_net", "cov_net_pu")),
      "Four-model registry mismatch"
    )
    pair <- c(
      cov = .59,
      network = round(4.16 * sqrt(2) / (qt(.975, 3) * sqrt(q * 399 / 400)), 2)
    )
    for (key in names(specs)) {
      sp <- specs[[key]]
      want <- if (key %in% c("cov", "cov_pu")) c(cov = .83) else pair
      lscar_assert(
        isTRUE(all.equal(sp$tau_scale, want, tolerance = 1e-12)) &&
          isTRUE(all.equal(sp$tau2_default, want^2, tolerance = 1e-12)),
        paste(key, "wrong prior assignment at q", q)
      )
      lscar_assert(
        identical(sp$use_pu, key %in% c("cov_pu", "cov_net_pu")),
        paste(key, "wrong PU flag")
      )
    }
  }
  invisible(TRUE)
}

lscar_linear_prior_selftest <- function(root) {
  src <- file.path(
    root,
    "simulation_linear_SCAR_all18/reference/sources/Bayesian_linear.R"
  )
  expr <- parse(src, keep.source = FALSE)
  index <- list()
  assignment_name <- function(z) {
    if (
      is.call(z) &&
        is.symbol(z[[1L]]) &&
        as.character(z[[1L]]) %in% c("<-", "=") &&
        length(z) >= 3L &&
        is.symbol(z[[2L]])
    ) {
      as.character(z[[2L]])
    } else {
      ""
    }
  }
  walk <- function(z) {
    if (!is.call(z)) {
      return(invisible(NULL))
    }
    head <- as.character(z[[1L]])[1L]
    nm <- assignment_name(z)
    if (
      nzchar(nm) &&
        is.call(z[[3L]]) &&
        identical(z[[3L]][[1L]], as.name("function"))
    ) {
      if (!is.null(index[[nm]])) {
        stop("Duplicate script-level function: ", nm)
      }
      index[[nm]] <<- z
      return(invisible(NULL))
    }
    if (head == "{" && length(z) >= 2L) {
      for (i in seq.int(2L, length(z))) {
        walk(z[[i]])
      }
    }
    if (head == "if" && length(z) >= 3L) {
      for (i in seq.int(3L, length(z))) {
        walk(z[[i]])
      }
    }
    invisible(NULL)
  }
  for (z in expr) {
    walk(z)
  }
  e <- new.env(parent = parent.env(.GlobalEnv))
  for (n in c(
    "round_prior",
    "logit",
    "solve_eta_beta_from_q95",
    "half_t_scale_from_geometry",
    "geometry_half_t_scales",
    "make_model_specs",
    "validate_groups"
  )) {
    if (is.null(index[[n]])) {
      stop("Missing source function: ", n)
    }
    eval(index[[n]], e)
  }
  settings <- c(
    "ROUND_ACTIVE_PRIORS",
    "PRIOR_ROUND_DIGITS",
    "TOTAL_SD_CENTER",
    "COV_VARIANCE_SHARE",
    "PAIR_VARIANCE_SHARES",
    "HALFT_DF",
    "MEAN_RMS_95_MULT",
    "BETA0_PREV_RANGE_95",
    "ETA_ALPHA",
    "ETA_Q95",
    "N_NODES",
    "NU_TAU",
    "beta0_low",
    "beta0_high",
    "beta0_prior_for_sampler",
    "eta_beta_unrounded",
    "eta_prior_for_sampler",
    "COVARIATE_MEAN_DIMENSION",
    "B_COV_FIXED",
    "HALF_T_Q95_UNIT",
    "EFFECTIVE_SD_CENTER_SINGLE",
    "EFFECTIVE_SD_CENTER_PAIR",
    "TAU_SCALE_COV_SINGLE_UNROUNDED",
    "TAU_SCALE_COV_SINGLE_ACTIVE",
    "TAU2_COV_SINGLE_DEFAULT",
    "TAU_SCALE_RULE"
  )
  nm <- vapply(expr, assignment_name, character(1))
  lscar_assert(all(settings %in% nm), "Missing source prior setting")
  for (k in which(nm %in% settings)) {
    eval(expr[[k]], e)
  }
  lscar_validate_linear_environment(e)
  bad <- lscar_linear_priors()
  bad$tau_scale_rule$tau_scale_cov_single <- .01
  lscar_assert(
    inherits(
      try(lscar_validate_linear_priors(bad, "negative test"), silent = TRUE),
      "try-error"
    ),
    "Prior mismatch not rejected"
  )
  cat("PASS: prior calibration of the Bayesian linear models.\n")
  invisible(TRUE)
}

lscar_linear_archive_checks <- function(archive) {
  d <- file.path(archive, "Bayesian_linear")
  cp <- list.files(
    file.path(d, "production_checkpoints"),
    pattern = "[.]rds$",
    full.names = TRUE
  )
  pp <- list.files(
    file.path(d, "pilot_checkpoints"),
    pattern = "[.]rds$",
    full.names = TRUE
  )
  pg <- list.files(
    file.path(d, "pilot_geometries"),
    pattern = "[.]rds$",
    full.names = TRUE
  )
  lscar_assert(
    length(cp) == 100L &&
      setequal(basename(cp), sprintf("replicate_%03d.rds", 1:100)),
    "Bayesian production inventory mismatch"
  )
  wanted <- as.vector(outer(
    sprintf("pilot_%02d_", 1:5),
    c("cov", "cov_pu", "cov_net", "cov_net_pu"),
    paste0
  ))
  lscar_assert(
    length(pp) == 20L && setequal(basename(pp), paste0(wanted, ".rds")),
    "Bayesian pilot inventory mismatch"
  )
  lscar_assert(
    length(pg) == 5L &&
      setequal(basename(pg), sprintf("pilot_%02d_geometry.rds", 1:5)),
    "Bayesian pilot geometries missing"
  )
  cf <- file.path(d, "global_pilot_calibration_5_datasets.rds")
  cal <- readRDS(cf)
  cfg <- readRDS(file.path(d, "Bayesian_linear_run_config.rds"))
  hash <- unname(tools::md5sum(cf))
  lscar_validate_linear_priors(
    cal$settings_signature$priors,
    "Global calibration"
  )
  lscar_assert(
    isTRUE(all.equal(
      cfg$prior_rule,
      lscar_linear_priors()$tau_scale_rule,
      tolerance = 1e-12
    )),
    "Run config prior rule mismatch"
  )
  lscar_assert(
    identical(as.character(cfg$global_calibration_hash), hash),
    "Run config references another calibration"
  )
  for (p in pp) {
    x <- readRDS(p)
    lscar_validate_linear_priors(x$settings_signature$priors, p)
    lscar_assert(
      identical(x$settings_signature, cal$settings_signature),
      paste("Pilot/global signature mismatch", p)
    )
    lscar_assert(
      x$settings_signature$pilot_burnin == 40000L &&
        x$settings_signature$pilot_tail_length == 5000L,
      paste("Pilot schedule mismatch", p)
    )
    lscar_assert(
      is.list(x$tail) && is.matrix(x$tail$tau2) && nrow(x$tail$tau2) == 5000L,
      paste("Incomplete pilot tail", p)
    )
  }
  for (p in cp) {
    x <- readRDS(p)
    lscar_validate_linear_priors(x$settings_signature$priors, p)
    lscar_assert(
      isTRUE(x$complete) &&
        identical(
          as.character(x$settings_signature$global_calibration_hash),
          hash
        ),
      paste("Production calibration mismatch", p)
    )
    fs <- x[["method_failures", exact = TRUE]]
    empty <- if (is.null(fs)) {
      TRUE
    } else if (is.data.frame(fs)) {
      nrow(fs) == 0L
    } else if (is.list(fs) && !is.object(fs) && is.null(dim(fs))) {
      length(fs) == 0L
    } else {
      FALSE
    }
    lscar_assert(empty, paste("Recorded or malformed failure", p))
  }
  cat(
    "PASS: Bayesian linear models: 100 production checkpoints, 20 pilots, five geometries and the prior calibration.\n"
  )
  invisible(TRUE)
}
