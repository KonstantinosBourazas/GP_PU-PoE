## Saving and reuse of completed fits.

original_model_signature <- function(model_spec, prior_bundle) {
  list(
    script_version = SCRIPT_VERSION,
    toy_data_hash = TOY_DATA_HASH,
    toy_seed = TOY_SEED,
    modularity_seed = MODULARITY_SEED,
    model_id = model_spec$model_id,
    method = model_spec$method,
    active_experts = model_spec$active_experts,
    use_pu = model_spec$use_pu,
    mcmc_seed = model_spec$mcmc_seed,
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
      covariate_background_mean = COVARIATE_BACKGROUND_MEAN,
      covariate_signal_mean = COVARIATE_SIGNAL_MEAN,
      common_sd = COMMON_SD,
      p_00 = P_00,
      p_01 = P_01,
      p_11 = P_11,
      modularity_q = modularity$q,
      modularity_selected_eigenvalues = modularity$selected_eigenvalues
    ),
    priors = list(
      beta0 = beta0_prior_for_sampler,
      eta = if (isTRUE(model_spec$use_pu)) {
        eta_prior_for_sampler
      } else {
        NULL
      },
      active_hyperpriors = prior_bundle$active_hyperpriors,
      tau_scale = prior_bundle$tau_scale,
      effective_sd_center = prior_bundle$effective_sd_center
    ),
    requested_interval = list(
      type = "HPD",
      package = "coda",
      function_name = "HPDinterval",
      probability = HPD_PROBABILITY
    )
  )
}

make_model_signature <- function(model_spec, prior_bundle) {
  list(
    original = original_model_signature(model_spec, prior_bundle),
    module_version = MODULE_VERSION,
    execution_mode = RUN_MODE,
    configuration = CONFIGURATION_SNAPSHOT,
    source_files_md5 = SOURCE_FILE_MD5,
    environment = ENVIRONMENT_FINGERPRINT
  )
}

atomic_save_rds <- function(object, file) {
  tmp <- tempfile(pattern = ".writing_", tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  ## Keep the previous complete file until the replacement has been serialized.
  if (file.exists(file) && !file.remove(file)) {
    stop("Cannot replace checkpoint: ", file)
  }
  if (!file.rename(tmp, file)) {
    stop("Cannot commit checkpoint: ", file)
  }
  invisible(file)
}

run_or_load_model <- function(model_spec) {
  prior_bundle <- model_prior_bundles[[model_spec$model_id]]
  model_signature <- make_model_signature(model_spec, prior_bundle)
  checkpoint_file <- file.path(
    MODEL_CHECKPOINT_DIR,
    paste0(model_spec$model_id, "_complete.rds")
  )
  if (
    isTRUE(RESUME_COMPLETED_MODELS) &&
      !isTRUE(FORCE_FRESH_RUN) &&
      file.exists(checkpoint_file)
  ) {
    candidate <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
    if (
      !is.list(candidate) ||
        !isTRUE(candidate$complete) ||
        !identical(candidate$model_signature, model_signature)
    ) {
      stop(
        "Checkpoint does not match the current code/settings/environment: ",
        checkpoint_file,
        ". Use a new --out folder or --fresh."
      )
    }
    d <- candidate$fit$draws
    expected <- c(N_RETAINED, N_NODES)
    if (
      !is.matrix(d$latent_probability) ||
        !is.matrix(d$observed_probability) ||
        !is.matrix(d$conditional_T_probability) ||
        !identical(dim(d$latent_probability), as.integer(expected)) ||
        !identical(dim(d$latent_probability), dim(d$observed_probability)) ||
        !identical(
          dim(d$latent_probability),
          dim(d$conditional_T_probability)
        ) ||
        any(!is.finite(d$latent_probability)) ||
        any(!is.finite(d$observed_probability)) ||
        any(!is.finite(d$conditional_T_probability))
    ) {
      stop("Checkpoint contains invalid probability draws: ", checkpoint_file)
    }
    cat("Reusing completed checkpoint for", model_spec$method, "\n")
    return(candidate)
  }
  set.seed(as.integer(model_spec$mcmc_seed))
  fit <- run_exact_toy_model_chain(
    model_spec = model_spec,
    prior_bundle = prior_bundle,
    Y = Y,
    X_cov = X_cov,
    A_network = A_network,
    Z_network = Z_network,
    network_geometry = network_geometry,
    m_beta = beta0_prior_for_sampler$mean,
    s_beta = beta0_prior_for_sampler$sd,
    a_eta = eta_prior_for_sampler$a_eta,
    b_eta = eta_prior_for_sampler$b_eta,
    n_blocked_burnin = BLOCKED_BURNIN,
    n_joint_burnin = JOINT_BURNIN,
    n_sample = POSTERIOR_ITERATIONS,
    thin = THIN,
    verbose = TRUE
  )
  summaries <- summarize_completed_model(fit, model_spec, prior_bundle)
  result <- list(
    complete = TRUE,
    model_signature = model_signature,
    model_spec = model_spec,
    prior_bundle = prior_bundle,
    fit = fit,
    node_summary = summaries$node_summary,
    eta_summary = summaries$eta_summary,
    diagnostics = summaries$diagnostics,
    prior_calibration = summaries$prior_calibration,
    completed_at = timestamp_now()
  )
  atomic_save_rds(result, checkpoint_file)
  result
}
