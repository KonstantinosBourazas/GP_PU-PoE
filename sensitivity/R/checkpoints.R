## Saving and reuse of the sensitivity fits.

sensitivity_signature <- function(fit_spec, model_spec, prior_bundle) {
  list(
    original = make_sensitivity_signature(fit_spec, model_spec, prior_bundle),
    module_version = SENSITIVITY_MODULE_VERSION,
    execution_mode = RUN_MODE,
    configuration = SENSITIVITY_CONFIGURATION,
    source_files_md5 = SENSITIVITY_SOURCE_MD5,
    environment = ENVIRONMENT_FINGERPRINT
  )
}

sensitivity_save_rds <- function(object, path) {
  ## Serialize and read-check the temporary file before replacing a summary.
  ## Complete MCMC checkpoints are not overwritten by the fitting workflow.
  tmp <- tempfile(".writing_sensitivity_", tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE, expand = FALSE), add = TRUE)
  saveRDS(object, tmp)
  invisible(readRDS(tmp))
  if (file.exists(path) && !file.remove(path)) {
    stop("Cannot replace output: ", path)
  }
  if (!file.rename(tmp, path)) {
    stop("Cannot commit output: ", path)
  }
  invisible(path)
}

validate_sensitivity_checkpoint <- function(x, signature = NULL) {
  if (
    !is.list(x) ||
      !isTRUE(x$complete) ||
      !is.list(x$fit) ||
      (!is.null(signature) && !identical(x$fit_signature, signature))
  ) {
    stop(
      "Incomplete checkpoint or code/settings/environment mismatch. ",
      "Use a new output folder."
    )
  }
  d <- x$fit$draws
  for (key in c(
    "latent_probability",
    "observed_probability",
    "conditional_T_probability"
  )) {
    z <- d[[key]]
    if (
      !is.matrix(z) ||
        !identical(dim(z), as.integer(c(N_RETAINED, N_NODES))) ||
        any(!is.finite(z)) ||
        any(z < 0 | z > 1)
    ) {
      stop("Invalid probability draws in checkpoint: ", key)
    }
  }
  if (
    length(d$eta) != N_RETAINED ||
      length(d$beta0) != N_RETAINED ||
      any(!is.finite(d$eta)) ||
      any(d$eta <= 0 | d$eta >= 1) ||
      any(!is.finite(d$beta0)) ||
      !identical(dim(d$tau), as.integer(c(N_RETAINED, 2L))) ||
      !identical(dim(d$theta), as.integer(c(N_RETAINED, 4L))) ||
      any(!is.finite(d$tau)) ||
      any(d$tau <= 0) ||
      any(!is.finite(d$theta)) ||
      any(d$theta <= 0)
  ) {
    stop("Invalid parameter draws in checkpoint.")
  }
  expected_phases <- c(
    blocked_warmup = BLOCKED_BURNIN,
    joint_warmup = JOINT_BURNIN,
    posterior_iterations = POSTERIOR_ITERATIONS,
    thinning = THIN,
    retained = N_RETAINED
  )
  if (!identical(as.integer(x$fit$phases), as.integer(expected_phases))) {
    stop("Checkpoint MCMC schedule does not match the selected mode.")
  }
  invisible(TRUE)
}

run_sensitivity_fit <- function(fit_id, load_only = FALSE) {
  fit_spec <- SENSITIVITY_SPECS[[fit_id]]
  if (is.null(fit_spec)) {
    stop("Unknown sensitivity specification: ", fit_id)
  }
  prior_bundle <- sensitivity_prior_bundles[[fit_id]]
  model_spec <- make_sensitivity_model_spec(fit_spec)
  signature <- sensitivity_signature(fit_spec, model_spec, prior_bundle)
  path <- file.path(FIT_CHECKPOINT_DIR, paste0(fit_id, "_complete.rds"))
  if (file.exists(path)) {
    x <- readRDS(path)
    validate_sensitivity_checkpoint(x, signature)
    cat("Reusing completed sensitivity fit", fit_id, "\n")
    return(x)
  }
  if (isTRUE(load_only)) {
    stop("Required completed fit not found: ", path)
  }
  cat("\n================ SENSITIVITY FIT:", fit_id, "================\n")
  cat(
    "GP SD:",
    fit_spec$total_sd_center,
    "; correlation targets:",
    paste(fit_spec$corr_pct_range, collapse = ", "),
    "; mean RMS multiplier:",
    fit_spec$mean_rms_95_mult,
    "\n"
  )
  cat(
    "beta0: N(",
    prior_bundle$beta0_prior$mean,
    ", ",
    prior_bundle$beta0_prior$sd,
    "^2); eta: Beta(",
    prior_bundle$eta_prior$a_eta,
    ", ",
    prior_bundle$eta_prior$b_eta,
    "); seed ",
    fit_spec$mcmc_seed,
    "\n",
    sep = ""
  )
  set.seed(as.integer(fit_spec$mcmc_seed))
  fit <- run_exact_toy_model_chain(
    model_spec = model_spec,
    prior_bundle = prior_bundle,
    Y = Y,
    X_cov = X_cov,
    A_network = A_network,
    Z_network = Z_network,
    network_geometry = network_geometry,
    m_beta = prior_bundle$beta0_prior$mean,
    s_beta = prior_bundle$beta0_prior$sd,
    a_eta = prior_bundle$eta_prior$a_eta,
    b_eta = prior_bundle$eta_prior$b_eta,
    n_blocked_burnin = BLOCKED_BURNIN,
    n_joint_burnin = JOINT_BURNIN,
    n_sample = POSTERIOR_ITERATIONS,
    thin = THIN,
    max_log_sigma_cov = log(fit_spec$max_sigma_cov),
    max_mean_sd_network = fit_spec$max_mean_sd_network,
    verbose = TRUE
  )
  summaries <- summarize_completed_model(fit, model_spec, prior_bundle)
  x <- list(
    complete = TRUE,
    fit_signature = signature,
    fit_spec = fit_spec,
    model_spec = model_spec,
    prior_bundle = prior_bundle,
    fit = fit,
    node_summary = add_fit_metadata(summaries$node_summary, fit_spec),
    eta_summary = add_fit_metadata(summaries$eta_summary, fit_spec),
    diagnostics = add_fit_metadata(summaries$diagnostics, fit_spec),
    ## The calibration table already holds the fit metadata columns.
    prior_calibration = prior_bundle$summary,
    completed_at = timestamp_now()
  )
  validate_sensitivity_checkpoint(x, signature)
  sensitivity_save_rds(x, path)
  x
}

compare_sensitivity_default <- function(default_result) {
  ## Compares the DEFAULT fit with the GP-PU-PoE fit of the illustrative example,
  ## when that fit exists. DEFAULT is always fitted here.
  path <- file.path(
    PROJECT_ROOT,
    "results",
    "illustration",
    RUN_MODE,
    "model_checkpoints",
    "GP_PU_POE_complete.rds"
  )
  report_path <- file.path(OUT_DIR, "default_vs_illustration.csv")
  if (!file.exists(path)) {
    write.csv(
      data.frame(
        check = "reference_checkpoint",
        status = "NOT_AVAILABLE",
        detail = "No fit of the illustrative example to compare with."
      ),
      report_path,
      row.names = FALSE
    )
    cat("DEFAULT: no fit of the illustrative example to compare with.\n")
    return(invisible(FALSE))
  }
  old <- readRDS(path)
  if (
    !is.list(old) ||
      !isTRUE(old$complete) ||
      is.null(old$model_signature$original)
  ) {
    stop("Not a completed fit of the illustrative example: ", path)
  }
  orig <- old$model_signature$original
  p <- default_result$prior_bundle
  prior_same <- identical(
    orig$priors$active_hyperpriors,
    p$active_hyperpriors
  ) &&
    identical(orig$priors$tau_scale, p$tau_scale) &&
    identical(orig$priors$effective_sd_center, p$effective_sd_center) &&
    identical(
      orig$priors$beta0,
      list(mean = p$beta0_prior$mean, sd = p$beta0_prior$sd)
    ) &&
    identical(
      orig$priors$eta,
      list(a_eta = p$eta_prior$a_eta, b_eta = p$eta_prior$b_eta)
    )
  keys <- setdiff(
    intersect(
      names(old$model_signature$configuration),
      names(CONFIGURATION_SNAPSHOT)
    ),
    c("SCRIPT_VERSION", "MODEL_SPECS")
  )
  settings_same <- identical(
    old$model_signature$configuration[keys],
    CONFIGURATION_SNAPSHOT[keys]
  )
  data_same <- identical(orig$toy_data_hash, TOY_DATA_HASH) &&
    identical(orig$data$modularity_q, modularity$q) &&
    identical(
      orig$data$modularity_selected_eigenvalues,
      modularity$selected_eigenvalues
    )
  env_same <- identical(
    old$model_signature$environment,
    ENVIRONMENT_FINGERPRINT
  )
  seed_same <- identical(
    as.integer(orig$mcmc_seed),
    as.integer(COMMON_MCMC_SEED)
  )
  fields <- c("draws", "state", "proposal", "acceptance", "phases")
  equal_fields <- vapply(
    fields,
    function(k) identical(default_result$fit[[k]], old$fit[[k]]),
    logical(1)
  )
  checks <- c(
    prior_bundle_values = prior_same,
    computational_configuration = settings_same,
    data_and_features = data_same,
    environment = env_same,
    MCMC_seed = seed_same,
    setNames(equal_fields, paste0("exact_", fields))
  )
  table <- data.frame(
    check = names(checks),
    identical = unname(checks),
    stringsAsFactors = FALSE
  )
  write.csv(table, report_path, row.names = FALSE)
  if (all(checks)) {
    cat("DEFAULT matches the fit of the illustrative example exactly.\n")
  } else if (
    prior_same && settings_same && data_same && env_same && seed_same
  ) {
    stop(
      "DEFAULT differs from the fit of the illustrative example although inputs, settings and environment match. ",
      "See default_vs_illustration.csv."
    )
  } else {
    warning(
      "The fit of the illustrative example has other inputs, settings or environment, so DEFAULT is not compared exactly. ",
      "See default_vs_illustration.csv.",
      call. = FALSE
    )
  }
  invisible(all(checks))
}
