## Settings of the prior sensitivity analysis. The data and MCMC settings come
## from config/illustration.R, so quick mode is that of the illustrative example.
SCRIPT_VERSION <- "toy_seed53_PU_PoE_prior_sensitivity_HPD80_v1_0"
COMMON_MCMC_SEED <- 530401L
TOP_K_VALUES <- c(10L, 20L, 30L)

SENSITIVITY_SPECS <- list(
  DEFAULT = make_sensitivity_spec(
    fit_id = "DEFAULT",
    varying_block = "default",
    prior_role = "default"
  ),
  SIGMA_CONSERVATIVE = make_sensitivity_spec(
    fit_id = "SIGMA_CONSERVATIVE",
    varying_block = "sigma",
    prior_role = "conservative",
    total_sd_center = 1
  ),
  SIGMA_DIFFUSE = make_sensitivity_spec(
    fit_id = "SIGMA_DIFFUSE",
    varying_block = "sigma",
    prior_role = "diffuse",
    total_sd_center = 4,
    max_sigma_cov = 8,
    max_mean_sd_network = 8
  ),
  ELL_CONSERVATIVE = make_sensitivity_spec(
    fit_id = "ELL_CONSERVATIVE",
    varying_block = "ell",
    prior_role = "conservative",
    corr_pct_range = c(0.05, 0.25)
  ),
  ELL_DIFFUSE = make_sensitivity_spec(
    fit_id = "ELL_DIFFUSE",
    varying_block = "ell",
    prior_role = "diffuse",
    corr_pct_range = c(0.10, 0.90)
  ),
  GAMMA_CONSERVATIVE = make_sensitivity_spec(
    fit_id = "GAMMA_CONSERVATIVE",
    varying_block = "gamma",
    prior_role = "conservative",
    mean_rms_95_mult = 2.08
  ),
  GAMMA_DIFFUSE = make_sensitivity_spec(
    fit_id = "GAMMA_DIFFUSE",
    varying_block = "gamma",
    prior_role = "diffuse",
    mean_rms_95_mult = 8.32
  ),
  BETA0_CONSERVATIVE = make_sensitivity_spec(
    fit_id = "BETA0_CONSERVATIVE",
    varying_block = "beta0",
    prior_role = "conservative",
    beta0_prev_range = c(0.005, 0.05)
  ),
  BETA0_DIFFUSE = make_sensitivity_spec(
    fit_id = "BETA0_DIFFUSE",
    varying_block = "beta0",
    prior_role = "diffuse",
    beta0_prev_range = c(0.02, 0.20)
  ),
  ETA_CONSERVATIVE = make_sensitivity_spec(
    fit_id = "ETA_CONSERVATIVE",
    varying_block = "eta",
    prior_role = "conservative",
    eta_prior_kind = "fixed_beta",
    eta_alpha = 1,
    eta_beta = 5
  ),
  ETA_DIFFUSE = make_sensitivity_spec(
    fit_id = "ETA_DIFFUSE",
    varying_block = "eta",
    prior_role = "diffuse",
    eta_prior_kind = "fixed_beta",
    eta_alpha = 1,
    eta_beta = 1
  )
)


JOINT_STRESS_SPECS <- list(
  JOINT_CONSERVATIVE = make_sensitivity_spec(
    fit_id = "JOINT_CONSERVATIVE",
    varying_block = "joint",
    prior_role = "conservative",
    total_sd_center = 1,
    corr_pct_range = c(0.05, 0.25),
    mean_rms_95_mult = 2.08,
    beta0_prev_range = c(0.005, 0.05),
    eta_prior_kind = "fixed_beta",
    eta_alpha = 1,
    eta_beta = 5
  ),
  JOINT_DIFFUSE = make_sensitivity_spec(
    fit_id = "JOINT_DIFFUSE",
    varying_block = "joint",
    prior_role = "diffuse",
    total_sd_center = 4,
    corr_pct_range = c(0.10, 0.90),
    mean_rms_95_mult = 8.32,
    beta0_prev_range = c(0.02, 0.20),
    eta_prior_kind = "fixed_beta",
    eta_alpha = 1,
    eta_beta = 1,
    max_sigma_cov = 8,
    max_mean_sd_network = 8
  )
)


SENSITIVITY_SPECS <- c(SENSITIVITY_SPECS, JOINT_STRESS_SPECS)
SENSITIVITY_RUN_ORDER <- names(SENSITIVITY_SPECS)

BASE_PUPOE_MODEL_SPEC <- list(
  method = "GP-PU-PoE",
  active_experts = c("cov", "network"),
  use_pu = TRUE,
  mcmc_seed = COMMON_MCMC_SEED
)


stopifnot(
  length(SENSITIVITY_SPECS) == 13L,
  !anyDuplicated(SENSITIVITY_RUN_ORDER),
  identical(SENSITIVITY_RUN_ORDER[1L], "DEFAULT")
)
