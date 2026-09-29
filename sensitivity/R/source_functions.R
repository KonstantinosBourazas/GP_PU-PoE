## Functions of the original sensitivity script, unchanged. Their line numbers
## in that script are in sensitivity/reference/SOURCE_FUNCTION_MAP.csv.

make_sensitivity_spec <- function(
    fit_id,
    varying_block,
    prior_role,
    total_sd_center = 2,
    corr_pct_range = c(0.10, 0.50),
    mean_rms_95_mult = 4.16,
    beta0_prev_range = c(0.01, 0.10),
    eta_prior_kind = "calibrated_q95",
    eta_alpha = 2,
    eta_beta = NA_real_,
    eta_q95 = 0.50,
    max_sigma_cov = 5,
    max_mean_sd_network = 5
) {
  list(
    fit_id = fit_id,
    varying_block = varying_block,
    prior_role = prior_role,
    total_sd_center = as.numeric(total_sd_center),
    corr_pct_range = as.numeric(corr_pct_range),
    mean_rms_95_mult = as.numeric(mean_rms_95_mult),
    beta0_prev_range = as.numeric(beta0_prev_range),
    eta_prior_kind = eta_prior_kind,
    eta_alpha = as.numeric(eta_alpha),
    eta_beta = as.numeric(eta_beta),
    eta_q95 = as.numeric(eta_q95),
    max_sigma_cov = as.numeric(max_sigma_cov),
    max_mean_sd_network = as.numeric(max_mean_sd_network),
    mcmc_seed = COMMON_MCMC_SEED
  )
}

calibrate_one_expert_prior <- function(
    expert_name,
    metric_fun,
    ell_bounds,
    effective_sd_center,
    corr_pct_range = CORR_PCT_RANGE_ALL,
    sigma_range_factor = SIGMA_RANGE_FACTOR
) {
  corr_pct_range <- as.numeric(corr_pct_range)
  sigma_range_factor <- as.numeric(sigma_range_factor)

  if (
    length(corr_pct_range) != 2L ||
    any(!is.finite(corr_pct_range)) ||
    corr_pct_range[1L] < 0 ||
    corr_pct_range[2L] > 1 ||
    corr_pct_range[1L] >= corr_pct_range[2L]
  ) {
    stop("corr_pct_range must be an increasing pair inside [0,1].")
  }

  if (
    length(sigma_range_factor) != 1L ||
    !is.finite(sigma_range_factor) ||
    sigma_range_factor <= 1
  ) {
    stop("sigma_range_factor must be one finite value larger than one.")
  }

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
    corr_pct_range * attainable_span

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
    effective_sd_center / sigma_range_factor
  ) / unit_kernel_mean_sd

  sigma_raw_high <- (
    effective_sd_center * sigma_range_factor
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
    target_corr_pct_low = corr_pct_range[1L],
    target_corr_pct_high = corr_pct_range[2L],
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
    effective_sd_target_low = effective_sd_center / sigma_range_factor,
    effective_sd_target_high = effective_sd_center * sigma_range_factor,
    sigma_range_factor = sigma_range_factor,
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

make_beta0_prior_from_range <- function(prevalence_range) {
  prevalence_range <- as.numeric(prevalence_range)

  if (
    length(prevalence_range) != 2L ||
    any(!is.finite(prevalence_range)) ||
    prevalence_range[1L] <= 0 ||
    prevalence_range[2L] >= 1 ||
    prevalence_range[1L] >= prevalence_range[2L]
  ) {
    stop("The beta0 prevalence range is invalid.")
  }

  lower_logit <- logit(prevalence_range[1L])
  upper_logit <- logit(prevalence_range[2L])

  list(
    mean = round_prior(0.5 * (lower_logit + upper_logit)),
    sd = round_prior(
      (upper_logit - lower_logit) /
        (2 * qnorm(0.975))
    ),
    prevalence_range = prevalence_range
  )
}

make_eta_prior_from_spec <- function(fit_spec) {
  if (identical(fit_spec$eta_prior_kind, "calibrated_q95")) {
    beta_value <- solve_eta_beta_from_q95(
      alpha = fit_spec$eta_alpha,
      q95_target = fit_spec$eta_q95,
      prob = 0.95
    )

    return(list(
      a_eta = as.numeric(fit_spec$eta_alpha),
      b_eta = round_prior(beta_value),
      construction = "calibrated_q95",
      q95_target = as.numeric(fit_spec$eta_q95)
    ))
  }

  if (identical(fit_spec$eta_prior_kind, "fixed_beta")) {
    if (
      !is.finite(fit_spec$eta_alpha) ||
      !is.finite(fit_spec$eta_beta) ||
      fit_spec$eta_alpha <= 0 ||
      fit_spec$eta_beta <= 0
    ) {
      stop("The fixed Beta prior for eta is invalid.")
    }

    return(list(
      a_eta = as.numeric(fit_spec$eta_alpha),
      b_eta = as.numeric(fit_spec$eta_beta),
      construction = "fixed_beta",
      q95_target = NA_real_
    ))
  }

  stop("Unknown eta_prior_kind: ", fit_spec$eta_prior_kind)
}

calibrate_sensitivity_fit_priors <- function(
    fit_spec,
    X_cov,
    Z_network,
    network_geometry
) {
  active_experts <- c("cov", "network")
  n_experts <- length(active_experts)

  effective_sd_center <- setNames(
    rep(fit_spec$total_sd_center / sqrt(n_experts), n_experts),
    active_experts
  )

  ## The structured-mean scale is a separate prior block. It is therefore
  ## referenced to the default expert SD budget even in the sigma-sensitivity
  ## fits, so changing TOTAL_SD_CENTER changes only the GP amplitude priors.
  structured_mean_sd_reference <- setNames(
    rep(TOTAL_SD_CENTER / sqrt(n_experts), n_experts),
    active_experts
  )

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
    } else {
      metric_fun <- network_metric_fun_factory(network_geometry)
      ell_bounds <- ELL_BOUNDS_NETWORK
      mean_design <- Z_network
    }

    calibration <- calibrate_one_expert_prior(
      expert_name = expert_name,
      metric_fun = metric_fun,
      ell_bounds = ell_bounds,
      effective_sd_center = effective_sd_center[expert_name],
      corr_pct_range = fit_spec$corr_pct_range,
      sigma_range_factor = SIGMA_RANGE_FACTOR
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
      fit_spec$mean_rms_95_mult * structured_mean_sd_reference[expert_name]
    ) / (
      half_t_q95_unit * b_factor
    )

    row <- calibration$row
    row$fit_id <- fit_spec$fit_id
    row$varying_block <- fit_spec$varying_block
    row$prior_role <- fit_spec$prior_role
    row$n_active_experts <- n_experts
    row$n_mean_features <- ncol(mean_design)
    row$Z_rms_factor <- b_factor
    row$mean_rms_95_mult <- fit_spec$mean_rms_95_mult
    row$total_sd_center <- fit_spec$total_sd_center
    row$structured_mean_sd_reference <-
      structured_mean_sd_reference[expert_name]
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
  beta0_prior <- make_beta0_prior_from_range(
    fit_spec$beta0_prev_range
  )
  eta_prior <- make_eta_prior_from_spec(fit_spec)

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
  summary_table$beta0_prior_mean <- beta0_prior$mean
  summary_table$beta0_prior_sd <- beta0_prior$sd
  summary_table$beta0_prev_low <- beta0_prior$prevalence_range[1L]
  summary_table$beta0_prev_high <- beta0_prior$prevalence_range[2L]
  summary_table$eta_prior_alpha <- eta_prior$a_eta
  summary_table$eta_prior_beta <- eta_prior$b_eta
  summary_table$max_sigma_cov <- fit_spec$max_sigma_cov
  summary_table$max_mean_sd_network <- fit_spec$max_mean_sd_network
  summary_table$active_prior_rounding <- ROUND_ACTIVE_PRIORS
  summary_table$active_prior_digits <- PRIOR_ROUND_DIGITS

  list(
    fit_id = fit_spec$fit_id,
    fit_spec = fit_spec,
    active_experts = active_experts,
    effective_sd_center = effective_sd_center,
    structured_mean_sd_reference = structured_mean_sd_reference,
    active_hyperpriors = active_hyperpriors,
    hyperpriors_unrounded = hyperpriors_unrounded,
    tau_scale = tau_scale,
    tau_scale_unrounded = tau_scale_unrounded,
    tau2_default = tau2_default,
    beta0_prior = beta0_prior,
    eta_prior = eta_prior,
    summary = summary_table,
    curves = curves
  )
}

make_sensitivity_model_spec <- function(fit_spec) {
  list(
    model_id = fit_spec$fit_id,
    method = "GP-PU-PoE",
    active_experts = BASE_PUPOE_MODEL_SPEC$active_experts,
    use_pu = TRUE,
    mcmc_seed = as.integer(fit_spec$mcmc_seed)
  )
}

make_sensitivity_signature <- function(
    fit_spec,
    model_spec,
    prior_bundle
) {
  list(
    script_version = SCRIPT_VERSION,
    toy_data_hash = TOY_DATA_HASH,
    toy_seed = TOY_SEED,
    modularity_seed = MODULARITY_SEED,
    fit_spec = fit_spec,
    model_spec = model_spec,
    MCMC = list(
      blocked_burnin = BLOCKED_BURNIN,
      joint_burnin = JOINT_BURNIN,
      posterior_iterations = POSTERIOR_ITERATIONS,
      thin = THIN,
      retained = N_RETAINED
    ),
    data = list(
      n_nodes = N_NODES,
      n_covariates = N_COVARIATES,
      true_zero_idx = TRUE_ZERO_IDX,
      hidden_positive_idx = HIDDEN_POS_IDX,
      observed_positive_idx = OBSERVED_POS_IDX,
      p_00 = P_00,
      p_01 = P_01,
      p_11 = P_11,
      modularity_q = modularity$q,
      modularity_selected_eigenvalues =
        modularity$selected_eigenvalues
    ),
    priors = list(
      beta0 = prior_bundle$beta0_prior,
      eta = prior_bundle$eta_prior,
      active_hyperpriors = prior_bundle$active_hyperpriors,
      tau_scale = prior_bundle$tau_scale,
      effective_sd_center = prior_bundle$effective_sd_center,
      max_sigma_cov = fit_spec$max_sigma_cov,
      max_mean_sd_network = fit_spec$max_mean_sd_network
    ),
    requested_interval = list(
      type = "HPD",
      package = "coda",
      function_name = "HPDinterval",
      probability = HPD_PROBABILITY
    )
  )
}

add_fit_metadata <- function(df, fit_spec) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  data.frame(
    fit_id = fit_spec$fit_id,
    varying_block = fit_spec$varying_block,
    prior_role = fit_spec$prior_role,
    df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

ordered_top_nodes <- function(probability, k) {
  node_index <- seq_along(probability)
  node_index[
    order(-probability, node_index)[seq_len(min(k, length(probability)))]
  ]
}

summarize_parameter_draws <- function(result_now) {
  fit_id_now <- result_now$fit_spec$fit_id
  spec_now <- result_now$fit_spec
  draws_now <- result_now$fit$draws

  parameter_vectors <- list(
    beta0 = as.numeric(draws_now$beta0),
    eta = as.numeric(draws_now$eta),
    tau_cov = as.numeric(draws_now$tau[, "tau_cov"]),
    tau_network = as.numeric(draws_now$tau[, "tau_network"]),
    ell_cov = as.numeric(draws_now$theta[, "ell_cov"]),
    sigma_cov = as.numeric(draws_now$theta[, "sigma_cov"]),
    ell_network = as.numeric(draws_now$theta[, "ell_network"]),
    sigma_network = as.numeric(draws_now$theta[, "sigma_network"])
  )

  rows <- lapply(
    names(parameter_vectors),
    function(parameter_name) {
      values <- parameter_vectors[[parameter_name]]
      interval <- hpd_interval_coda(
        values,
        probability = HPD_PROBABILITY
      )

      data.frame(
        fit_id = fit_id_now,
        varying_block = spec_now$varying_block,
        prior_role = spec_now$prior_role,
        parameter = parameter_name,
        n_draws = length(values),
        posterior_mean = mean(values),
        posterior_sd = stats::sd(values),
        HPD80_lower = interval["lower"],
        HPD80_upper = interval["upper"],
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(rbind, rows)
}

effective_size_matrix <- function(draw_matrix) {
  draw_matrix <- as.matrix(draw_matrix)

  if (
    nrow(draw_matrix) < 2L ||
    ncol(draw_matrix) < 1L ||
    any(!is.finite(draw_matrix))
  ) {
    stop("ESS received an invalid posterior-draw matrix.")
  }

  as.numeric(
    coda::effectiveSize(
      coda::mcmc(draw_matrix)
    )
  )
}
