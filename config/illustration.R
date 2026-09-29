## Settings of the illustrative example: the full settings of the paper, then
## the shorter quick mode. RUN_MODE is set by the workflow.

SCRIPT_VERSION <- "toy_seed53_four_GP_models_HPD80_v1_0"

## Settings of the synthetic dataset.
TOY_SEED <- 53L
N_NODES <- 100L
N_COVARIATES <- 10L

TRUE_ZERO_IDX <- 1:70
HIDDEN_POS_IDX <- 71:80
OBSERVED_POS_IDX <- 81:100
TRUE_POS_IDX <- 71:100

COVARIATE_BACKGROUND_MEAN <- 0.00
COVARIATE_SIGNAL_MEAN <- 0.50
COMMON_SD <- 1.00

P_00 <- 0.10
P_01 <- 0.05
P_11 <- 0.16

## A separate seed is used by RSpectra after the data have been generated, so
## it cannot change X, A, T or Y.
MODULARITY_SEED <- 53043L
MOD_EIG_TOL <- 1e-8
MOD_MAX_Q <- Inf
MOD_K_MAX <- 20L

## MCMC lengths of the paper.
BLOCKED_BURNIN <- 10000L
JOINT_BURNIN <- 20000L
POSTERIOR_ITERATIONS <- 150000L
THIN <- 30L
N_RETAINED <- POSTERIOR_ITERATIONS %/% THIN
HPD_PROBABILITY <- 0.80


## Prior construction, as in the full PU-PoE code.
ROUND_ACTIVE_PRIORS <- TRUE
PRIOR_ROUND_DIGITS <- 2L

CORR_PCT_RANGE_ALL <- c(0.10, 0.50)
ELL_INTERVAL_PROBS <- c(0.025, 0.975)
SIGMA_INTERVAL_PROBS <- c(0.025, 0.975)
N_ELL_GRID_COARSE <- 60L
ROOT_TOL <- 1e-8

TOTAL_SD_CENTER <- 2
SIGMA_RANGE_FACTOR <- 2

HALFT_DF <- 3
MEAN_RMS_95_MULT <- 4.16
NU_TAU <- HALFT_DF
SAMPLE_TAU_GROUPS <- TRUE

BETA0_PREV_RANGE_95 <- c(0.01, 0.10)
ETA_ALPHA <- 2
ETA_Q95 <- 0.50

## Proposal and adaptation settings of the sampler.
KAPPA_INIT <- 1
ADAPT_BLOCK_CORRELATIONS <- TRUE
COR_ADAPT_START <- 2000L
COR_ADAPT_INTERVAL <- 250L
COR_ADAPT_WINDOW <- 5000L
COR_SHRINKAGE <- 0.05
COR_MAX_ABS <- 0.95

EMP_COV_TAIL <- 5000L
JOINT_COR_SHRINKAGE <- 0.10
JOINT_COR_MAX_ABS <- 0.95
F_REFRESH_AFTER_JOINT <- TRUE

ADAPT <- TRUE
ADAPT_START <- 2000L
ADAPT_INTERVAL <- 25L
TARGET_F <- 0.55
TARGET_THETA <- 0.25
TARGET_JOINT <- 0.25
ADAPT_RATE_DELTA <- 0.05
ADAPT_RATE_THETA <- 0.05
ADAPT_RATE_JOINT_THETA <- 0.5
T0_ADAPT <- 10
POWER_ADAPT <- 0.60
DELTA_INIT <- 3.5
DELTA_MIN <- 1e-4
DELTA_MAX <- 50
KAPPA_MIN <- 1e-8
KAPPA_MAX <- 100
JOINT_THETA_SCALE_MIN <- 1e-4
JOINT_THETA_SCALE_MAX <- 20
COV_JITTER <- 1e-6

ETA_SLICE_WIDTH <- 0.50
INITIAL_ETA <- 0.20
ETA_EPS <- 1e-10
JITTER <- 1e-8

MAX_LOG_SIGMA_COV <- log(5)
MAX_MEAN_SD_NETWORK <- 5
MAX_LOG_SIGMA_NETWORK <- log(1e3)
MIN_LOG_SIGMA <- log(1e-3)
LOGL_COV_BOUNDS <- c(log(1e-3), log(1e3))
LOGL_NETWORK_BOUNDS <- c(log(1e-4), log(1e2))
ELL_BOUNDS_COV <- exp(LOGL_COV_BOUNDS)
ELL_BOUNDS_NETWORK <- exp(LOGL_NETWORK_BOUNDS)

PI_SPIKE <- 0
SPIKE_LOG_SIG_MEAN <- -6
SPIKE_LOG_SIG_SD <- 0.5

BLOCKED_PROGRESS_EVERY <- 1000L
JOINT_PROGRESS_EVERY <- 2000L
SAMPLE_PROGRESS_EVERY <- 5000L

MODEL_SPECS <- list(
  list(
    model_id = "GP_COV",
    method = "GP-Cov",
    active_experts = "cov",
    use_pu = FALSE,
    mcmc_seed = 530101L
  ),
  list(
    model_id = "GP_NET",
    method = "GP-Net",
    active_experts = "network",
    use_pu = FALSE,
    mcmc_seed = 530201L
  ),
  list(
    model_id = "GP_POE",
    method = "GP-PoE",
    active_experts = c("cov", "network"),
    use_pu = FALSE,
    mcmc_seed = 530301L
  ),
  list(
    model_id = "GP_PU_POE",
    method = "GP-PU-PoE",
    active_experts = c("cov", "network"),
    use_pu = TRUE,
    mcmc_seed = 530401L
  )
)

## Quick mode changes only the lengths of the runs. Data, priors and seeds stay
## the same, and the draws only check that the code runs.
if (identical(RUN_MODE, "quick")) {
  BLOCKED_BURNIN <- 200L
  JOINT_BURNIN <- 300L
  POSTERIOR_ITERATIONS <- 600L
  THIN <- 3L
  ADAPT_START <- 50L
  ADAPT_INTERVAL <- 25L
  COR_ADAPT_START <- 100L
  COR_ADAPT_INTERVAL <- 50L
  COR_ADAPT_WINDOW <- 200L
  EMP_COV_TAIL <- 200L
  BLOCKED_PROGRESS_EVERY <- 100L
  JOINT_PROGRESS_EVERY <- 100L
  SAMPLE_PROGRESS_EVERY <- 200L
} else if (identical(RUN_MODE, "custom")) {
  ## A new run of the full design with another MCMC seed or sampling length.
  ## Data, priors, warm up and adaptation stay those of the full run. With
  ## CUSTOM_MCMC_SEED = 53 the model seeds are those of the paper.
  POSTERIOR_ITERATIONS <- as.integer(CUSTOM_ITERATIONS)
  THIN <- as.integer(CUSTOM_THIN)
  for (k in seq_along(MODEL_SPECS)) {
    MODEL_SPECS[[k]]$mcmc_seed <- as.integer(CUSTOM_MCMC_SEED) *
      10000L +
      100L * k +
      1L
  }
} else if (!identical(RUN_MODE, "full")) {
  stop("RUN_MODE must be 'quick', 'full' or 'custom'.")
}

N_RETAINED <- POSTERIOR_ITERATIONS %/% THIN

validate_illustration_schedule(
  BLOCKED_BURNIN,
  JOINT_BURNIN,
  POSTERIOR_ITERATIONS,
  THIN,
  paper_schedule = identical(RUN_MODE, "full")
)

if (ADAPT_START + ADAPT_INTERVAL - 1L > min(BLOCKED_BURNIN, JOINT_BURNIN)) {
  stop("The adaptation schedule would not update in both warm-up phases.")
}

if (COR_ADAPT_START > BLOCKED_BURNIN || EMP_COV_TAIL < 2L) {
  stop("The proposal-correlation schedule is invalid.")
}

CONFIGURATION_KEYS <- c(
  "SCRIPT_VERSION",
  "TOY_SEED",
  "N_NODES",
  "N_COVARIATES",
  "TRUE_ZERO_IDX",
  "HIDDEN_POS_IDX",
  "OBSERVED_POS_IDX",
  "TRUE_POS_IDX",
  "COVARIATE_BACKGROUND_MEAN",
  "COVARIATE_SIGNAL_MEAN",
  "COMMON_SD",
  "P_00",
  "P_01",
  "P_11",
  "MODULARITY_SEED",
  "MOD_EIG_TOL",
  "MOD_MAX_Q",
  "MOD_K_MAX",
  "BLOCKED_BURNIN",
  "JOINT_BURNIN",
  "POSTERIOR_ITERATIONS",
  "THIN",
  "N_RETAINED",
  "HPD_PROBABILITY",
  "ROUND_ACTIVE_PRIORS",
  "PRIOR_ROUND_DIGITS",
  "CORR_PCT_RANGE_ALL",
  "ELL_INTERVAL_PROBS",
  "SIGMA_INTERVAL_PROBS",
  "N_ELL_GRID_COARSE",
  "ROOT_TOL",
  "TOTAL_SD_CENTER",
  "SIGMA_RANGE_FACTOR",
  "HALFT_DF",
  "MEAN_RMS_95_MULT",
  "NU_TAU",
  "SAMPLE_TAU_GROUPS",
  "BETA0_PREV_RANGE_95",
  "ETA_ALPHA",
  "ETA_Q95",
  "KAPPA_INIT",
  "ADAPT_BLOCK_CORRELATIONS",
  "COR_ADAPT_START",
  "COR_ADAPT_INTERVAL",
  "COR_ADAPT_WINDOW",
  "COR_SHRINKAGE",
  "COR_MAX_ABS",
  "EMP_COV_TAIL",
  "JOINT_COR_SHRINKAGE",
  "JOINT_COR_MAX_ABS",
  "F_REFRESH_AFTER_JOINT",
  "ADAPT",
  "ADAPT_START",
  "ADAPT_INTERVAL",
  "TARGET_F",
  "TARGET_THETA",
  "TARGET_JOINT",
  "ADAPT_RATE_DELTA",
  "ADAPT_RATE_THETA",
  "ADAPT_RATE_JOINT_THETA",
  "T0_ADAPT",
  "POWER_ADAPT",
  "DELTA_INIT",
  "DELTA_MIN",
  "DELTA_MAX",
  "KAPPA_MIN",
  "KAPPA_MAX",
  "JOINT_THETA_SCALE_MIN",
  "JOINT_THETA_SCALE_MAX",
  "COV_JITTER",
  "ETA_SLICE_WIDTH",
  "INITIAL_ETA",
  "ETA_EPS",
  "JITTER",
  "MAX_LOG_SIGMA_COV",
  "MAX_MEAN_SD_NETWORK",
  "MAX_LOG_SIGMA_NETWORK",
  "MIN_LOG_SIGMA",
  "LOGL_COV_BOUNDS",
  "LOGL_NETWORK_BOUNDS",
  "ELL_BOUNDS_COV",
  "ELL_BOUNDS_NETWORK",
  "PI_SPIKE",
  "SPIKE_LOG_SIG_MEAN",
  "SPIKE_LOG_SIG_SD",
  "BLOCKED_PROGRESS_EVERY",
  "JOINT_PROGRESS_EVERY",
  "SAMPLE_PROGRESS_EVERY",
  "MODEL_SPECS"
)
