## Extracted from GP_four_model_illustrtation.R. See docs/SOURCE_MAP.csv.
## Original mathematical operations and their order are retained.

solve_eta_beta_from_q95 <- function(
    alpha,
    q95_target,
    prob = 0.95
) {
  if (
    !is.finite(alpha) || alpha <= 0 ||
    !is.finite(q95_target) ||
    q95_target <= 0 ||
    q95_target >= 1
  ) {
    stop("Invalid eta-prior target.")
  }

  objective <- function(beta) {
    qbeta(
      prob,
      shape1 = alpha,
      shape2 = beta
    ) - q95_target
  }

  lo <- 1e-6
  hi <- 1

  while (objective(hi) > 0) {
    hi <- hi * 2
    if (hi > 1e6) {
      stop("Could not bracket the eta-prior beta parameter.")
    }
  }

  uniroot(objective, interval = c(lo, hi))$root
}

fit_lognormal_interval <- function(
    q_low,
    q_high,
    probs = c(0.025, 0.975)
) {
  q_low <- as.numeric(q_low)
  q_high <- as.numeric(q_high)

  if (
    length(q_low) != 1L ||
    length(q_high) != 1L ||
    !is.finite(q_low) ||
    !is.finite(q_high) ||
    q_low <= 0 ||
    q_high <= 0
  ) {
    stop("q_low and q_high must be positive finite scalars.")
  }

  qlo <- min(q_low, q_high)
  qhi <- max(q_low, q_high)
  p_low <- probs[1L]
  p_high <- probs[2L]

  if (
    !is.finite(p_low) ||
    !is.finite(p_high) ||
    p_low <= 0 ||
    p_high >= 1 ||
    p_low >= p_high
  ) {
    stop("Invalid lognormal interval probabilities.")
  }

  tau <- (
    log(qhi) - log(qlo)
  ) / (
    qnorm(p_high) - qnorm(p_low)
  )

  mu <- log(qlo) - qnorm(p_low) * tau

  list(
    mu = mu,
    tau = tau,
    q_low = qlo,
    q_high = qhi,
    p_low = p_low,
    p_high = p_high,
    median = exp(mu),
    mean = exp(mu + 0.5 * tau^2)
  )
}

make_covariate_metric_fun <- function(X_std) {
  X_std <- as.matrix(X_std)
  pair_d2 <- as.numeric(stats::dist(X_std))^2
  expected_pairs <- nrow(X_std) * (nrow(X_std) - 1) / 2

  if (length(pair_d2) != expected_pairs) {
    stop("Unexpected covariate-pair count in prior calibration.")
  }

  function(ell) {
    ell <- as.numeric(ell)

    if (
      length(ell) != 1L ||
      !is.finite(ell) ||
      ell <= 0
    ) {
      return(list(
        mean_corr = NA_real_,
        mean_sd = NA_real_
      ))
    }

    list(
      mean_corr = mean(exp(-pair_d2 / (2 * ell^2))),
      mean_sd = 1
    )
  }
}

network_metric_fun_factory <- function(network_geometry) {
  n <- nrow(network_geometry$A)
  total_pairs <- n * (n - 1) / 2
  components <- network_geometry$components

  function(ell) {
    ell <- as.numeric(ell)

    if (
      length(ell) != 1L ||
      !is.finite(ell) ||
      ell <= 0
    ) {
      return(list(
        mean_corr = NA_real_,
        mean_sd = NA_real_
      ))
    }

    corr_sum <- 0
    sd_sum <- 0

    for (component in components) {
      s <- component$size

      if (component$type == "isolate") {
        diag_value <- exp(-1 / (2 * ell^2))
        sd_sum <- sd_sum + sqrt(diag_value)
        next
      }

      if (component$type == "clique") {
        q <- exp(-component$lambda_orth / (2 * ell^2))
        diag_value <- q * (1 - 1 / s) + 1 / s
        off_value <- (1 - q) / s
        corr_value <- off_value / diag_value

        corr_sum <- corr_sum + choose(s, 2) * corr_value
        sd_sum <- sd_sum + s * sqrt(diag_value)
        next
      }

      gamma <- exp(-component$lambda / (2 * ell^2))
      diag_value <- as.vector(component$U2 %*% gamma)
      diag_value <- pmax(diag_value, 1e-12)
      sd_vector <- sqrt(diag_value)
      inverse_sd <- 1 / sd_vector
      U_trans_inverse_sd <- as.vector(
        crossprod(component$U, inverse_sd)
      )
      full_corr_sum <- sum(gamma * U_trans_inverse_sd^2)
      off_corr_sum <- 0.5 * (full_corr_sum - s)

      if (
        is.finite(off_corr_sum) &&
        off_corr_sum < 0 &&
        abs(off_corr_sum) < 1e-8
      ) {
        off_corr_sum <- 0
      }

      corr_sum <- corr_sum + off_corr_sum
      sd_sum <- sd_sum + sum(sd_vector)
    }

    mean_corr <- corr_sum / total_pairs
    mean_corr <- min(max(mean_corr, 0), 1)

    list(
      mean_corr = mean_corr,
      mean_sd = sd_sum / n
    )
  }
}

build_corr_curve <- function(
    metric_fun,
    ell_bounds,
    n_grid = N_ELL_GRID_COARSE
) {
  ell_grid <- exp(seq(
    log(ell_bounds[1L]),
    log(ell_bounds[2L]),
    length.out = as.integer(n_grid)
  ))

  mean_corr <- numeric(length(ell_grid))
  mean_sd <- numeric(length(ell_grid))

  for (ii in seq_along(ell_grid)) {
    value <- metric_fun(ell_grid[ii])
    mean_corr[ii] <- value$mean_corr
    mean_sd[ii] <- value$mean_sd
  }

  ok <- is.finite(ell_grid) &
    is.finite(mean_corr) &
    is.finite(mean_sd)

  data.frame(
    ell = ell_grid[ok],
    mean_corr = mean_corr[ok],
    mean_sd = mean_sd[ok],
    stringsAsFactors = FALSE
  )
}

invert_corr_curve_refined <- function(
    curve,
    metric_fun,
    target_corr,
    expert_name
) {
  curve <- curve[
    is.finite(curve$ell) &
      is.finite(curve$mean_corr),
    ,
    drop = FALSE
  ]

  curve <- curve[order(curve$ell), , drop = FALSE]

  corr_range <- range(curve$mean_corr, na.rm = TRUE)
  target_used <- min(
    max(target_corr, corr_range[1L]),
    corr_range[2L]
  )

  difference <- curve$mean_corr - target_used
  exact_idx <- which(abs(difference) <= 1e-12)

  if (length(exact_idx) > 0L) {
    ell_hat <- curve$ell[exact_idx[1L]]
  } else {
    crossing_idx <- which(
      difference[-nrow(curve)] *
        difference[-1L] <= 0
    )

    root_result <- NULL

    if (length(crossing_idx) > 0L) {
      crossing_score <- abs(difference[crossing_idx]) +
        abs(difference[crossing_idx + 1L])
      kk <- crossing_idx[which.min(crossing_score)]
      log_interval <- log(c(
        curve$ell[kk],
        curve$ell[kk + 1L]
      ))

      root_result <- tryCatch(
        uniroot(
          function(log_ell) {
            metric_fun(exp(log_ell))$mean_corr - target_used
          },
          interval = log_interval,
          tol = ROOT_TOL
        ),
        error = function(e) NULL
      )
    }

    if (!is.null(root_result)) {
      ell_hat <- exp(root_result$root)
    } else {
      ordering <- order(curve$mean_corr, curve$ell)
      corr_ordered <- curve$mean_corr[ordering]
      ell_ordered <- curve$ell[ordering]
      keep <- !duplicated(round(corr_ordered, 12))
      corr_unique <- corr_ordered[keep]
      ell_unique <- ell_ordered[keep]

      if (length(corr_unique) >= 2L) {
        ell_hat <- as.numeric(approx(
          x = corr_unique,
          y = ell_unique,
          xout = target_used,
          rule = 2,
          ties = "ordered"
        )$y)
      } else {
        ell_hat <- curve$ell[which.min(abs(difference))]
      }
    }
  }

  if (!is.finite(ell_hat) || ell_hat <= 0) {
    stop(
      "Failed to invert the correlation curve for ",
      expert_name,
      "."
    )
  }

  list(
    ell = ell_hat,
    target_requested = target_corr,
    target_used = target_used,
    achieved_corr = metric_fun(ell_hat)$mean_corr,
    attainable_min = corr_range[1L],
    attainable_max = corr_range[2L]
  )
}

calibrate_one_expert_prior <- function(
    expert_name,
    metric_fun,
    ell_bounds,
    effective_sd_center
) {
  curve <- build_corr_curve(
    metric_fun = metric_fun,
    ell_bounds = ell_bounds,
    n_grid = N_ELL_GRID_COARSE
  )

  attainable <- range(curve$mean_corr, na.rm = TRUE)
  attainable_span <- diff(attainable)

  if (!is.finite(attainable_span) || attainable_span <= 1e-10) {
    stop(
      "The attainable correlation range is degenerate for ",
      expert_name,
      "."
    )
  }

  targets <- attainable[1L] +
    CORR_PCT_RANGE_ALL * attainable_span

  inversion_1 <- invert_corr_curve_refined(
    curve = curve,
    metric_fun = metric_fun,
    target_corr = targets[1L],
    expert_name = expert_name
  )

  inversion_2 <- invert_corr_curve_refined(
    curve = curve,
    metric_fun = metric_fun,
    target_corr = targets[2L],
    expert_name = expert_name
  )

  ell_low <- min(inversion_1$ell, inversion_2$ell)
  ell_high <- max(inversion_1$ell, inversion_2$ell)

  prior_ell <- fit_lognormal_interval(
    q_low = ell_low,
    q_high = ell_high,
    probs = ELL_INTERVAL_PROBS
  )

  ell_median <- exp(prior_ell$mu)
  unit_kernel_mean_sd <- metric_fun(ell_median)$mean_sd

  if (!is.finite(unit_kernel_mean_sd) || unit_kernel_mean_sd <= 0) {
    stop(
      "Invalid unit-amplitude mean marginal SD for ",
      expert_name,
      "."
    )
  }

  sigma_raw_low <- (
    effective_sd_center / SIGMA_RANGE_FACTOR
  ) / unit_kernel_mean_sd

  sigma_raw_high <- (
    effective_sd_center * SIGMA_RANGE_FACTOR
  ) / unit_kernel_mean_sd

  prior_sigma <- fit_lognormal_interval(
    q_low = sigma_raw_low,
    q_high = sigma_raw_high,
    probs = SIGMA_INTERVAL_PROBS
  )

  hyperprior_unrounded <- list(
    log_ell_mean = prior_ell$mu,
    log_ell_sd = prior_ell$tau,
    log_sig_mean = prior_sigma$mu,
    log_sig_sd = prior_sigma$tau
  )

  row <- data.frame(
    expert = expert_name,
    attainable_corr_min = attainable[1L],
    attainable_corr_max = attainable[2L],
    attainable_corr_span = attainable_span,
    target_corr_pct_low = CORR_PCT_RANGE_ALL[1L],
    target_corr_pct_high = CORR_PCT_RANGE_ALL[2L],
    target_corr_low = targets[1L],
    target_corr_high = targets[2L],
    achieved_corr_at_first_target = inversion_1$achieved_corr,
    achieved_corr_at_second_target = inversion_2$achieved_corr,
    ell_interval_low = ell_low,
    ell_interval_high = ell_high,
    log_ell_mean_unrounded = prior_ell$mu,
    log_ell_sd_unrounded = prior_ell$tau,
    ell_median = prior_ell$median,
    unit_kernel_mean_marginal_sd = unit_kernel_mean_sd,
    effective_sd_target_center = effective_sd_center,
    effective_sd_target_low = effective_sd_center / SIGMA_RANGE_FACTOR,
    effective_sd_target_high = effective_sd_center * SIGMA_RANGE_FACTOR,
    sigma_raw_interval_low = sigma_raw_low,
    sigma_raw_interval_high = sigma_raw_high,
    log_sig_mean_unrounded = prior_sigma$mu,
    log_sig_sd_unrounded = prior_sigma$tau,
    sigma_raw_median = prior_sigma$median,
    effective_sigma_median_check =
      prior_sigma$median * unit_kernel_mean_sd,
    stringsAsFactors = FALSE
  )

  list(
    hyperprior_unrounded = hyperprior_unrounded,
    row = row,
    curve = curve
  )
}

activate_hyperprior <- function(hyperprior_unrounded) {
  required_fields <- c(
    "log_ell_mean",
    "log_ell_sd",
    "log_sig_mean",
    "log_sig_sd"
  )

  if (!all(required_fields %in% names(hyperprior_unrounded))) {
    stop("A geometry-specific hyperprior is incomplete.")
  }

  output <- lapply(
    hyperprior_unrounded[required_fields],
    as.numeric
  )

  if (isTRUE(ROUND_ACTIVE_PRIORS)) {
    output <- lapply(output, round_prior)
  }

  if (output$log_ell_sd <= 0) {
    output$log_ell_sd <- hyperprior_unrounded$log_ell_sd
  }

  if (output$log_sig_sd <= 0) {
    output$log_sig_sd <- hyperprior_unrounded$log_sig_sd
  }

  if (
    any(!is.finite(unlist(output))) ||
    output$log_ell_sd <= 0 ||
    output$log_sig_sd <= 0
  ) {
    stop("The active geometry-specific hyperprior is invalid.")
  }

  output
}

calibrate_model_priors <- function(
    model_spec,
    X_cov,
    Z_network,
    network_geometry
) {
  active_experts <- model_spec$active_experts
  n_experts <- length(active_experts)

  if (n_experts < 1L || n_experts > 2L) {
    stop("The toy illustration supports one or two active experts.")
  }

  ## A one-expert ablation gets the full total SD budget. In a two-expert model
  ## the two experts split the total residual variance equally, as in the full
  ## model.
  effective_sd_center <- if (n_experts == 1L) {
    setNames(TOTAL_SD_CENTER, active_experts)
  } else {
    setNames(
      rep(TOTAL_SD_CENTER / sqrt(n_experts), n_experts),
      active_experts
    )
  }

  active_hyperpriors <- list()
  hyperpriors_unrounded <- list()
  summary_rows <- list()
  curves <- list()
  tau_scale_unrounded <- numeric(n_experts)
  names(tau_scale_unrounded) <- active_experts

  half_t_q95_unit <- qt(
    (0.95 + 1) / 2,
    df = HALFT_DF
  )

  for (expert_name in active_experts) {
    if (identical(expert_name, "cov")) {
      metric_fun <- make_covariate_metric_fun(X_cov)
      ell_bounds <- ELL_BOUNDS_COV
      mean_design <- X_cov
    } else if (identical(expert_name, "network")) {
      metric_fun <- network_metric_fun_factory(network_geometry)
      ell_bounds <- ELL_BOUNDS_NETWORK
      mean_design <- Z_network
    } else {
      stop("Unknown expert: ", expert_name)
    }

    calibration <- calibrate_one_expert_prior(
      expert_name = expert_name,
      metric_fun = metric_fun,
      ell_bounds = ell_bounds,
      effective_sd_center = effective_sd_center[expert_name]
    )

    active_hyperpriors[[expert_name]] <- activate_hyperprior(
      calibration$hyperprior_unrounded
    )

    hyperpriors_unrounded[[expert_name]] <-
      calibration$hyperprior_unrounded

    b_factor <- sqrt(
      sum(mean_design^2) /
        nrow(mean_design)
    )

    if (!is.finite(b_factor) || b_factor <= 0) {
      stop("The structured-mean RMS factor is invalid for ", expert_name, ".")
    }

    tau_scale_unrounded[expert_name] <- (
      MEAN_RMS_95_MULT * effective_sd_center[expert_name]
    ) / (
      half_t_q95_unit * b_factor
    )

    row <- calibration$row
    row$model_id <- model_spec$model_id
    row$method <- model_spec$method
    row$n_active_experts <- n_experts
    row$n_mean_features <- ncol(mean_design)
    row$Z_rms_factor <- b_factor
    row$tau_scale_unrounded <- tau_scale_unrounded[expert_name]
    summary_rows[[expert_name]] <- row
    curves[[expert_name]] <- calibration$curve
  }

  tau_scale <- if (isTRUE(ROUND_ACTIVE_PRIORS)) {
    round_prior(tau_scale_unrounded)
  } else {
    tau_scale_unrounded
  }
  names(tau_scale) <- active_experts

  if (any(!is.finite(tau_scale)) || any(tau_scale <= 0)) {
    stop("At least one active half-t scale is invalid.")
  }

  tau2_default <- tau_scale^2

  summary_table <- do.call(rbind, summary_rows)
  row.names(summary_table) <- NULL
  summary_table$tau_scale_active <- unname(
    tau_scale[summary_table$expert]
  )
  summary_table$tau2_default <- unname(
    tau2_default[summary_table$expert]
  )
  summary_table$log_ell_mean_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_ell_mean
    },
    numeric(1)
  )
  summary_table$log_ell_sd_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_ell_sd
    },
    numeric(1)
  )
  summary_table$log_sig_mean_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_sig_mean
    },
    numeric(1)
  )
  summary_table$log_sig_sd_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_sig_sd
    },
    numeric(1)
  )
  summary_table$active_prior_rounding <- ROUND_ACTIVE_PRIORS
  summary_table$active_prior_digits <- PRIOR_ROUND_DIGITS

  list(
    model_id = model_spec$model_id,
    method = model_spec$method,
    active_experts = active_experts,
    effective_sd_center = effective_sd_center,
    active_hyperpriors = active_hyperpriors,
    hyperpriors_unrounded = hyperpriors_unrounded,
    tau_scale = tau_scale,
    tau_scale_unrounded = tau_scale_unrounded,
    tau2_default = tau2_default,
    summary = summary_table,
    curves = curves
  )
}

