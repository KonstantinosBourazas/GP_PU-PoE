## Extracted from GP_four_model_illustrtation.R. See docs/SOURCE_MAP.csv.
## Original mathematical operations and their order are retained.

hpd_interval_coda <- function(x, probability = HPD_PROBABILITY) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    stop("At least two finite draws are required for an HPD interval.")
  }

  if (!is.finite(probability) || probability <= 0 || probability >= 1) {
    stop("The HPD probability must lie strictly between zero and one.")
  }

  ## coda::HPDinterval may receive a constant chain (e.g. Pr(T=1|Y=1)=1 in
  ## the PU model). In that special case the exact HPD interval is degenerate.
  if (max(x) - min(x) <= 10 * .Machine$double.eps) {
    value <- x[1L]
    return(c(lower = value, upper = value))
  }

  interval <- coda::HPDinterval(
    coda::as.mcmc(x),
    prob = probability
  )

  c(
    lower = as.numeric(interval[1L, "lower"]),
    upper = as.numeric(interval[1L, "upper"])
  )
}

hpd_matrix_coda <- function(
    draw_matrix,
    probability = HPD_PROBABILITY
) {
  draw_matrix <- as.matrix(draw_matrix)

  if (
    nrow(draw_matrix) < 2L ||
    ncol(draw_matrix) < 1L ||
    any(!is.finite(draw_matrix))
  ) {
    stop("The draw matrix supplied to coda HPD calculation is invalid.")
  }

  intervals <- t(vapply(
    seq_len(ncol(draw_matrix)),
    function(column_index) {
      hpd_interval_coda(
        draw_matrix[, column_index],
        probability = probability
      )
    },
    numeric(2)
  ))

  colnames(intervals) <- c("lower", "upper")
  intervals
}

roc_auc <- function(scores, labels) {
  scores <- as.numeric(scores)
  labels <- as.integer(labels)
  ok <- is.finite(scores) & labels %in% c(0L, 1L)
  scores <- scores[ok]
  labels <- labels[ok]

  n_positive <- sum(labels == 1L)
  n_negative <- sum(labels == 0L)

  if (n_positive == 0L || n_negative == 0L) {
    return(NA_real_)
  }

  ranks_ascending <- rank(scores, ties.method = "average")
  sum_ranks_positive <- sum(ranks_ascending[labels == 1L])

  (
    sum_ranks_positive -
      n_positive * (n_positive + 1) / 2
  ) / (
    n_positive * n_negative
  )
}

overall_logloss <- function(
    probabilities,
    truth,
    eps = 1e-8
) {
  probabilities <- pmin(
    pmax(as.numeric(probabilities), eps),
    1 - eps
  )
  truth <- as.numeric(truth)

  -mean(
    truth * log(probabilities) +
      (1 - truth) * log(1 - probabilities)
  )
}

overall_brier <- function(probabilities, truth) {
  mean(
    (as.numeric(probabilities) - as.numeric(truth))^2
  )
}

summarize_completed_model <- function(
    fit,
    model_spec,
    prior_bundle
) {
  latent_draws <- as.matrix(fit$draws$latent_probability)
  observed_draws <- as.matrix(fit$draws$observed_probability)
  conditional_T_draws <- as.matrix(
    fit$draws$conditional_T_probability
  )

  expected_dimension <- c(N_RETAINED, N_NODES)

  if (
    !all(dim(latent_draws) == expected_dimension) ||
    !all(dim(observed_draws) == expected_dimension) ||
    !all(dim(conditional_T_draws) == expected_dimension)
  ) {
    stop(
      model_spec$method,
      " did not return the configured probability-matrix dimensions."
    )
  }

  latent_mean <- colMeans(latent_draws)
  observed_mean <- colMeans(observed_draws)
  conditional_T_mean <- colMeans(conditional_T_draws)

  latent_hpd <- hpd_matrix_coda(
    latent_draws,
    probability = HPD_PROBABILITY
  )

  observed_hpd <- hpd_matrix_coda(
    observed_draws,
    probability = HPD_PROBABILITY
  )

  conditional_T_hpd <- hpd_matrix_coda(
    conditional_T_draws,
    probability = HPD_PROBABILITY
  )

  node_summary <- data.frame(
    model_id = model_spec$model_id,
    method = model_spec$method,
    active_experts = paste(
      model_spec$active_experts,
      collapse = " + "
    ),
    use_PU_likelihood = isTRUE(model_spec$use_pu),
    node = seq_len(N_NODES),
    T = T_TRUE,
    Y = Y,
    node_group = NODE_GROUP,

    ## Primary figure-ready posterior summaries: p_i = plogis(f_i).
    posterior_predictive_mean = latent_mean,
    HPD80_lower = latent_hpd[, "lower"],
    HPD80_upper = latent_hpd[, "upper"],
    HPD80_width = latent_hpd[, "upper"] - latent_hpd[, "lower"],

    latent_probability_mean = latent_mean,
    latent_HPD80_lower = latent_hpd[, "lower"],
    latent_HPD80_upper = latent_hpd[, "upper"],
    latent_HPD80_width = latent_hpd[, "upper"] - latent_hpd[, "lower"],

    ## Probability of the observed label under each fitted model.
    observed_probability_mean = observed_mean,
    observed_HPD80_lower = observed_hpd[, "lower"],
    observed_HPD80_upper = observed_hpd[, "upper"],
    observed_HPD80_width =
      observed_hpd[, "upper"] - observed_hpd[, "lower"],

    ## Draw-by-draw Pr(T_i=1 | Y_i, parameters). This differs from p_i only
    ## for GP-PU-PoE.
    conditional_T_probability_mean = conditional_T_mean,
    conditional_T_HPD80_lower = conditional_T_hpd[, "lower"],
    conditional_T_HPD80_upper = conditional_T_hpd[, "upper"],
    conditional_T_HPD80_width =
      conditional_T_hpd[, "upper"] - conditional_T_hpd[, "lower"],

    stringsAsFactors = FALSE
  )

  unlabeled_idx <- which(Y == 0L)
  hidden_truth_unlabeled <- as.integer(T_TRUE[unlabeled_idx] == 1L)

  eta_summary <- if (isTRUE(model_spec$use_pu)) {
    eta_hpd <- hpd_interval_coda(
      fit$draws$eta,
      probability = HPD_PROBABILITY
    )

    data.frame(
      model_id = model_spec$model_id,
      method = model_spec$method,
      eta_true = length(HIDDEN_POS_IDX) / length(TRUE_POS_IDX),
      eta_posterior_mean = mean(fit$draws$eta),
      eta_posterior_sd = stats::sd(fit$draws$eta),
      eta_HPD80_lower = eta_hpd["lower"],
      eta_HPD80_upper = eta_hpd["upper"],
      eta_HPD80_width = eta_hpd["upper"] - eta_hpd["lower"],
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      model_id = model_spec$model_id,
      method = model_spec$method,
      eta_true = length(HIDDEN_POS_IDX) / length(TRUE_POS_IDX),
      eta_posterior_mean = NA_real_,
      eta_posterior_sd = NA_real_,
      eta_HPD80_lower = NA_real_,
      eta_HPD80_upper = NA_real_,
      eta_HPD80_width = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  blocked_acceptance <- fit$acceptance$blocked

  diagnostics <- data.frame(
    model_id = model_spec$model_id,
    method = model_spec$method,
    active_experts = paste(
      model_spec$active_experts,
      collapse = " + "
    ),
    use_PU_likelihood = isTRUE(model_spec$use_pu),
    n_nodes = N_NODES,
    n_covariates = N_COVARIATES,
    n_network_features = ncol(Z_network),
    blocked_burnin = fit$phases["blocked_warmup"],
    joint_burnin = fit$phases["joint_warmup"],
    posterior_iterations = fit$phases["posterior_iterations"],
    thinning = fit$phases["thinning"],
    retained_draws = fit$phases["retained"],
    acc_blocked_f = unname(blocked_acceptance["f"]),
    acc_blocked_cov = if ("cov" %in% names(blocked_acceptance)) {
      unname(blocked_acceptance["cov"])
    } else {
      NA_real_
    },
    acc_blocked_network = if ("network" %in% names(blocked_acceptance)) {
      unname(blocked_acceptance["network"])
    } else {
      NA_real_
    },
    acc_joint_tune = fit$acceptance$joint_tune,
    acc_f_refresh_tune = fit$acceptance$f_refresh_tune,
    acc_joint_sampling = fit$acceptance$joint_sampling,
    acc_f_refresh_sampling = fit$acceptance$f_refresh_sampling,
    eta_move_rate_sampling = fit$acceptance$eta_moved_sampling,
    delta_final = fit$proposal$delta_final,
    joint_theta_scale_final = fit$proposal$joint_theta_scale_final,
    runtime_sec = fit$runtime_sec,
    overall_AUC_latent_mean = roc_auc(latent_mean, T_TRUE),
    hidden_AUC_latent_mean = roc_auc(
      latent_mean[unlabeled_idx],
      hidden_truth_unlabeled
    ),
    overall_LogLoss_latent_mean = overall_logloss(
      latent_mean,
      T_TRUE
    ),
    overall_Brier_latent_mean = overall_brier(
      latent_mean,
      T_TRUE
    ),
    mean_HPD80_width = mean(node_summary$HPD80_width),
    mean_HPD80_width_true_zero = mean(
      node_summary$HPD80_width[node_summary$node_group == "true zero"]
    ),
    mean_HPD80_width_hidden_positive = mean(
      node_summary$HPD80_width[
        node_summary$node_group == "hidden positive"
      ]
    ),
    mean_HPD80_width_observed_positive = mean(
      node_summary$HPD80_width[
        node_summary$node_group == "observed positive"
      ]
    ),
    stringsAsFactors = FALSE
  )

  list(
    node_summary = node_summary,
    eta_summary = eta_summary,
    diagnostics = diagnostics,
    prior_calibration = prior_bundle$summary
  )
}

get_complete_model_result <- function(model_id_now) {
  
  ## During the original full run, model_results is still in memory.
  if (
    exists("model_results", inherits = TRUE)
  ) {
    model_results_now <- get(
      "model_results",
      inherits = TRUE
    )
    
    candidate <- model_results_now[[model_id_now]]
    
    if (
      is.list(candidate) &&
      !is.null(candidate$fit) &&
      !is.null(candidate$fit$draws)
    ) {
      return(candidate)
    }
  }
  
  ## Otherwise read the complete model checkpoint.
  checkpoint_file <- file.path(
    MODEL_CHECKPOINT_DIR,
    paste0(
      model_id_now,
      "_complete.rds"
    )
  )
  
  if (!file.exists(checkpoint_file)) {
    stop(
      "Missing completed checkpoint for ",
      model_id_now,
      ": ",
      checkpoint_file
    )
  }
  
  candidate <- readRDS(
    checkpoint_file
  )
  
  if (
    !is.list(candidate) ||
    !isTRUE(candidate$complete) ||
    is.null(candidate$fit) ||
    is.null(candidate$fit$draws)
  ) {
    stop(
      "The checkpoint for ",
      model_id_now,
      " does not contain complete posterior draws."
    )
  }
  
  candidate
}

effective_size_matrix <- function(draw_matrix) {
  
  draw_matrix <- as.matrix(
    draw_matrix
  )
  
  if (
    nrow(draw_matrix) < 2L ||
    ncol(draw_matrix) < 1L ||
    any(!is.finite(draw_matrix))
  ) {
    stop(
      "ESS received an invalid posterior-draw matrix."
    )
  }
  
  as.numeric(
    coda::effectiveSize(
      coda::mcmc(
        draw_matrix
      )
    )
  )
}

