## Checks of the inputs before any fit.
check_sensitivity_inputs <- function() {
  expected_ids <- c(
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
  stopifnot(
    identical(SENSITIVITY_RUN_ORDER, expected_ids),
    identical(names(sensitivity_prior_bundles), expected_ids),
    N_NODES == 100L,
    ncol(X_cov) == 10L,
    ncol(Z_network) == 1L,
    sum(T_TRUE) == 30L,
    sum(Y) == 20L,
    sum(T_TRUE == 1L & Y == 0L) == 10L,
    sum(A_network) / 2 == 390,
    COMMON_MCMC_SEED == 530401L
  )
  ## The 13 prior bundles must equal those of the original script.
  original <- parse(
    file.path(
      PROJECT_ROOT,
      "sensitivity",
      "reference",
      "PU_PoE_prior_sensitivity_final.R"
    ),
    keep.source = FALSE
  )
  ref <- new.env(parent = as.environment("package:stats"))
  for (key in c(CONFIGURATION_KEYS, "COMMON_MCMC_SEED", "TOP_K_VALUES")) {
    assign(key, get(key, inherits = TRUE), envir = ref)
  }
  for (expr in original) {
    if (
      is.call(expr) &&
        identical(expr[[1L]], as.name("<-")) &&
        length(expr) == 3L
    ) {
      rhs <- expr[[3L]]
      if (is.call(rhs) && identical(rhs[[1L]], as.name("function"))) {
        eval(expr, ref)
      }
    }
  }
  for (key in c("SENSITIVITY_SPECS", "JOINT_STRESS_SPECS")) {
    ix <- which(vapply(
      original,
      function(x) {
        is.call(x) &&
          length(x) == 3L &&
          identical(x[[1L]], as.name("<-")) &&
          identical(x[[2L]], as.name(key))
      },
      logical(1)
    ))
    if (length(ix) != 1L) {
      stop("Cannot locate original specifications for ", key)
    }
    eval(original[[ix]], ref)
  }
  original_specs <- c(ref$SENSITIVITY_SPECS, ref$JOINT_STRESS_SPECS)
  if (!identical(SENSITIVITY_SPECS, original_specs)) {
    stop("The specifications differ from the original script.")
  }
  for (id in expected_ids) {
    bundle <- ref$calibrate_sensitivity_fit_priors(
      original_specs[[id]],
      X_cov,
      Z_network,
      network_geometry
    )
    if (!identical(bundle, sensitivity_prior_bundles[[id]])) {
      stop("Prior calibration differs from the original script for ", id)
    }
  }
  ## DEFAULT must have the priors of the illustrative example.
  core <- new.env(parent = as.environment("package:stats"))
  sys.source(
    file.path(PROJECT_ROOT, "R", "illustration", "workflow.R"),
    envir = core
  )
  ie <- core$new_illustration_environment(PROJECT_ROOT, mode = RUN_MODE)
  ms <- ie$MODEL_SPECS[[4L]]
  ip <- ie$calibrate_model_priors(ms, X_cov, Z_network, network_geometry)
  dp <- sensitivity_prior_bundles$DEFAULT
  for (k in c(
    "active_hyperpriors",
    "tau_scale",
    "tau2_default",
    "effective_sd_center"
  )) {
    if (!identical(ip[[k]], dp[[k]])) {
      stop("DEFAULT prior does not match illustration: ", k)
    }
  }
  stopifnot(
    dp$beta0_prior$mean == -3.40,
    dp$beta0_prior$sd == 0.61,
    dp$eta_prior$a_eta == 2,
    dp$eta_prior$b_eta == 6.39,
    SENSITIVITY_SPECS$SIGMA_DIFFUSE$max_sigma_cov == 8,
    SENSITIVITY_SPECS$SIGMA_DIFFUSE$max_mean_sd_network == 8,
    SENSITIVITY_SPECS$JOINT_DIFFUSE$max_sigma_cov == 8,
    SENSITIVITY_SPECS$JOINT_DIFFUSE$max_mean_sd_network == 8,
    identical(
      dp$tau_scale,
      sensitivity_prior_bundles$SIGMA_CONSERVATIVE$tau_scale
    ),
    identical(dp$tau_scale, sensitivity_prior_bundles$SIGMA_DIFFUSE$tau_scale)
  )
  cat(
    "Inputs checked: 13 specifications and priors, data, DEFAULT priors, diffuse support and mean blocks.\n"
  )
  invisible(TRUE)
}
