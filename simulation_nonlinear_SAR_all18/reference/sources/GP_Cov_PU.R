## ============================================================================
## SIMULATION STUDY: SELF-CALIBRATING GP-COV + PU MODEL UNDER SAR
## ----------------------------------------------------------------------------
## MODEL
## -----
## This is the covariate-only Gaussian-process model with positive-unlabeled
## observation noise:
##
##   f = beta0 + X gamma_cov + g_cov,
##   g_cov ~ GP(0, K_RBF),
##   Y_i | f_i, eta ~ Bernoulli((1 - eta) * plogis(f_i)).
##
## The latent class probability is p_i = plogis(f_i). In every pilot and
## production replication, 30 of the 45 true positives are observed and
## the remaining 15 are hidden among the 370 observations with Y = 0.
##
## COVARIATE-DRIVEN SAR LABEL MECHANISM
## ------------------------------------
## The latent truth and mixed-signal geometry are identical to the SCAR study:
## 355 true zeros and 45 true positives, split into 15 covariate-only,
## 15 network-only and 15 both-signal positives.
##
## In every pilot and production replication, exactly 15 of the 45 positives are
## hidden without replacement with weights proportional to
##
##   exp(-lambda * z_i),  lambda = 2,
##
## where z_i is the oracle covariate fraud-likeness score standardized within
## the 45 positives and clipped to [-3,3]. The oracle score is used only to
## generate missing labels and is never supplied to the fitted GP-Cov + PU model.
##
## SELF-CONTAINED PRIOR CALIBRATION
## --------------------------------
## This script has no dependency on the real-data application and no dependency
## on an external simulation-prior file.
##
## For every pilot and every production replication, before MCMC:
##
##   1. the current 400 x 10 standardized covariate geometry is used to
##      calibrate the prior for log(ell_cov);
##   2. the central 95% prior interval for ell_cov spans the values attaining
##      10% and 50% of the mean-correlation range attainable by that X;
##   3. the single GP expert receives effective marginal-SD center 2, with
##      central 95% range [1, 4], which calibrates log(sigma_cov);
##   4. gamma_cov | tau_cov ~ N(0, tau_cov^2 I_10), with
##      tau_cov ~ half-t_3(0, s_cov);
##   5. s_cov is fixed from the ten standardized covariate columns so that the
##      95th percentile of the expert-level RMS structured-mean contribution
##      equals 4.16 times the expert SD budget 2.
##
## The beta0 and eta priors are created once inside this script:
##
##   - beta0 is calibrated from the central 95% baseline-prevalence range
##     [0.01, 0.10];
##   - eta ~ Beta(2, b_eta), with its 95th percentile equal to 0.50.
##
## FIVE GLOBAL PILOTS
## ------------------
## Five independent SAR pilot data sets are used only for transferable MCMC
## initialization and proposal geometry. Since the ell prior changes with the
## current X geometry, theta and tau states are pooled on prior-relative scales
## and mapped back to the exact prior of each production replication. The eta
## state is pooled directly because its prior is fixed. The current-data latent
## f vector is never transferred between data sets.
##
## PRODUCTION AND EVALUATION
## -------------------------
##   - 100 independent SAR replications;
##   - exact dense 400 x 400 covariate GP;
##   - 15,000 production burn-in iterations;
##   - 40,000 posterior-sampling iterations, thinned every 20;
##   - all overall and hidden point metrics;
##   - 80% and 95% predictive-label coverage and interval scores.
##
## For every retained posterior draw, p_i^(s) = plogis(f_i^(s)). Point scores
## are posterior means of p_i. Predictive latent-label draws, used for the
## coverage and interval scores, are T_i ~ Bernoulli(p_i^(s)) for every unit
## and every retained draw. Predictive intervals use
## discrete empirical quantiles with type = 1 and are compared with the fixed
## latent truth. Hidden positives always have T_i = 1.
## ============================================================================

options(stringsAsFactors = FALSE)

## ============================================================================
## 0. PACKAGE CHECKS
## ============================================================================

required_namespaces <- c("Matrix")

missing_namespaces <- required_namespaces[
  !vapply(required_namespaces, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_namespaces) > 0L) {
  stop(
    "Install the required package(s) before running this script: ",
    paste(missing_namespaces, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Matrix)
})

## ============================================================================
## 1. USER SETTINGS
## ============================================================================

SCRIPT_VERSION <- "gp_cov_pu_self_calibrating_priors_SAR_lambda2_all_metrics_v1_0"
EVALUATION_VERSION <- "SAR_lambda2_node_predictive_labels_80_95_v2_0"

## Production batch to run in the current R session.
PRODUCTION_REPLICATES <- 1:25L
BATCH_ONLY <- TRUE

MASTER_SEED <- 20260718L
N_PILOT <- 5L
N_REP <- 100L

PILOT_BURNIN <- 40000L
PILOT_TAIL_LENGTH <- 5000L
PILOT_BLOCKED_FRAC <- 1 / 3

PRODUCTION_BURNIN <- 15000L
PRODUCTION_SAMPLING <- 40000L
PRODUCTION_THIN <- 20L
PRODUCTION_BLOCKED_FRAC <- 1 / 3

MIN_BLOCKED_WARMUP <- 5000L

## --------------------------------------------------------------------------
## Exact geometry-specific prior calibration
## --------------------------------------------------------------------------

## Exact unrounded values are always saved. The active MCMC prior is rounded to
## two decimals after calibration, matching the convention used in the other
## final simulation scripts.
ROUND_ACTIVE_PRIORS <- TRUE
PRIOR_ROUND_DIGITS <- 2L

CORR_PCT_RANGE_ALL <- c(0.10, 0.50)
ELL_INTERVAL_PROBS <- c(0.025, 0.975)
SIGMA_INTERVAL_PROBS <- c(0.025, 0.975)
N_ELL_GRID_COARSE <- 60L
ROOT_TOL <- 1e-8

## GP-Cov + PU contains one residual GP expert, so it receives the full SD budget 2.
TOTAL_SD_CENTER <- 2
EFFECTIVE_SD_CENTER_COV <- TOTAL_SD_CENTER
SIGMA_RANGE_FACTOR <- 2

## Structured-mean half-t calibration.
HALFT_DF <- 3
MEAN_RMS_95_MULT <- 4.16
NU_TAU <- HALFT_DF
SAMPLE_TAU <- TRUE

## Fixed beta0 and eta prior rules.
BETA0_PREV_RANGE_95 <- c(0.01, 0.10)
ETA_ALPHA <- 2
ETA_Q95 <- 0.50
ETA_SLICE_WIDTH <- 0.50
INITIAL_ETA <- 0.20
ETA_EPS <- 1e-10

## Hard bounds used by both calibration and MCMC.
MAX_LOG_SIGMA_COV <- log(5)
MIN_LOG_SIGMA <- log(1e-3)
LOGL_COV_BOUNDS <- c(log(1e-3), log(1e3))
ELL_BOUNDS_COV <- exp(LOGL_COV_BOUNDS)

## Proposal and adaptation controls.
KAPPA_COV_INIT <- 1
PROPOSAL_COR_INIT <- diag(2)

ADAPT_BLOCK_CORRELATION_PILOTS <- TRUE
ADAPT_BLOCK_CORRELATION_PRODUCTION <- FALSE
COR_ADAPT_START <- 2000L
COR_ADAPT_INTERVAL <- 250L
COR_ADAPT_WINDOW <- 5000L
COR_SHRINKAGE <- 0.05
COR_MAX_ABS <- 0.95
GLOBAL_COR_SHRINKAGE <- 0.10

JOINT_THETA_SCALE_INIT_PILOT <- 1
JOINT_THETA_SCALE_INIT_PRODUCTION <- 1
ADAPT_RATE_JOINT_THETA <- 0.5
F_REFRESH_AFTER_JOINT <- TRUE

ADAPT <- TRUE
ADAPT_START <- 2000L
ADAPT_INTERVAL <- 25L
TARGET_F <- 0.55
TARGET_THETA <- 0.25
TARGET_JOINT <- 0.25
ADAPT_RATE_DELTA <- 0.05
ADAPT_RATE_THETA <- 0.05
T0_ADAPT <- 10
POWER_ADAPT <- 0.60
DELTA_MIN <- 1e-4
DELTA_MAX <- 50
KAPPA_MIN <- 1e-8
KAPPA_MAX <- 100
JOINT_THETA_SCALE_MIN <- 1e-4
JOINT_THETA_SCALE_MAX <- 20
EMP_COV_TAIL <- 5000L
COV_JITTER <- 1e-6
JITTER <- 1e-8

PILOT_PROGRESS_EVERY <- 5000L
PRODUCTION_BURNIN_PROGRESS_EVERY <- 1000L
PRODUCTION_SAMPLE_PROGRESS_EVERY <- 5000L

HIDDEN_TOP_K <- 15L
OVERALL_RECALL_K <- 45L
PREDICTIVE_INTERVAL_LEVELS <- c(0.80, 0.95)
PREDICTIVE_QUANTILE_TYPE <- 1L
PROBABILITY_EPS <- 1e-8

OUT_DIR <- file.path(
  getwd(),
  "SIMULATION_100_GP_COV_PU_SELF_CALIBRATING_PRIORS_SAR_LAMBDA_2"
)

PILOT_CHECKPOINT_DIR <- file.path(OUT_DIR, "pilot_checkpoints")
PRODUCTION_CHECKPOINT_DIR <- file.path(OUT_DIR, "production_checkpoints")
PILOT_GEOMETRY_DIR <- file.path(OUT_DIR, "pilot_geometries")
PRIOR_CALIBRATION_DIR <- file.path(OUT_DIR, "prior_calibrations")

for (path_now in c(
  OUT_DIR,
  PILOT_CHECKPOINT_DIR,
  PRODUCTION_CHECKPOINT_DIR,
  PILOT_GEOMETRY_DIR,
  PRIOR_CALIBRATION_DIR
)) {
  dir.create(path_now, showWarnings = FALSE, recursive = TRUE)
}

GLOBAL_CALIBRATION_FILE <- file.path(
  OUT_DIR,
  "GP_Cov_PU_global_pilot_calibration_5_datasets.rds"
)

GLOBAL_CALIBRATION_RDATA <- file.path(
  OUT_DIR,
  "GP_Cov_PU_global_pilot_calibration_5_datasets.RData"
)

RESUME_PILOTS <- TRUE
RESUME_PRODUCTION <- TRUE
FORCE_FRESH_PILOTS <- FALSE
FORCE_FRESH_PRODUCTION <- FALSE
FAIL_IF_PRODUCTION_FAILS <- TRUE
SAVE_PILOT_GEOMETRIES <- TRUE

## Use the same X matrices as the competing-method simulation whenever
## they have already been saved. Otherwise regenerate them from the same seeds.
USE_SAVED_COMPETING_GEOMETRIES <- TRUE

COMPETING_OUT_DIR <- file.path(
  getwd(),
  "SIMULATION_100_COMPETING_METHODS_MIXED_SIGNAL_SAR_LAMBDA_2"
)

COMPETING_GEOMETRY_DIR <- file.path(
  COMPETING_OUT_DIR,
  "generated_geometries"
)

## Deterministic, independent MCMC and predictive-label streams.
PILOT_MCMC_SEED_BASE <- 710000000L
PRODUCTION_MCMC_SEED_BASE <- 810000000L
PRODUCTION_PREDICTIVE_SEED_BASE <- 910000000L

## ============================================================================
## 1A. FIXED LATENT TRUTH, MIXED-SIGNAL GEOMETRY AND SAR SETTINGS
## ============================================================================

N_NODES <- 400L
TRUE_ZERO_IDX <- 1:355
TRUE_POS_IDX <- 356:400

## These are the positive signal groups used by the SCAR design.
## Their hidden/observed status is regenerated in every replication.
COV_ONLY_IDX <- c(
  356:360,
  371:380
)

NET_ONLY_IDX <- c(
  361:365,
  381:390
)

BOTH_IDX <- c(
  366:370,
  391:400
)

COVARIATE_SIGNAL_IDX <- c(
  COV_ONLY_IDX,
  BOTH_IDX
)

NETWORK_SIGNAL_IDX <- c(
  NET_ONLY_IDX,
  BOTH_IDX
)

N_TRUE_POS <- length(TRUE_POS_IDX)
N_HIDDEN <- 15L
N_OBSERVED_POS <- N_TRUE_POS - N_HIDDEN
N_TRUE_ZERO <- length(TRUE_ZERO_IDX)
N_UNLABELED <- N_TRUE_ZERO + N_HIDDEN

TRUE_T <- integer(N_NODES)
TRUE_T[TRUE_POS_IDX] <- 1L

POSITIVE_SIGNAL_GROUP <- rep(
  "true zero",
  N_NODES
)

POSITIVE_SIGNAL_GROUP[COV_ONLY_IDX] <- "covariates only"
POSITIVE_SIGNAL_GROUP[NET_ONLY_IDX] <- "network only"
POSITIVE_SIGNAL_GROUP[BOTH_IDX] <- "both"

LAMBDA_SAR <- 2.00
SAR_Z_CLIP <- 3.00

SCENARIO_NAME <- paste0(
  "Fixed mixed-signal covariate-driven SAR design, lambda = ",
  format(LAMBDA_SAR, trim = TRUE)
)

METHOD_NAME <- "GP-Cov + PU"
MODEL_KEY <- "gp_cov_pu"

## ============================================================================
## 1B. COVARIATE DGP SETTINGS
## ============================================================================

MU_LINEAR <- 0.50
SD_LINEAR <- 1.00

BLOCK2_X3_OFFSET <- 2.60
BLOCK2_X4_CENTER <- 1.00
BLOCK2_X4_OFFSET <- 1.40
BLOCK2_ZERO_SD_X3 <- 2.00
BLOCK2_ZERO_SD_X4 <- 1.20
BLOCK2_POS_SD_X3 <- 0.90
BLOCK2_POS_SD_X4 <- 0.70

DISC_RADIUS <- 2.00
RING_RADIUS <- 1.50
RING_SD <- 0.25

PARABOLA_COEF <- 0.35
DELTA_PARABOLA <- 0.80
SD_PARABOLA_ZERO <- 0.40
SD_PARABOLA_POS <- 0.35

RHO_ZERO <- 0.85
RHO_POS <- -0.85

## ============================================================================

## ============================================================================
## 1C. REPRODUCIBLE SEED TABLES
## ============================================================================

RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

set.seed(MASTER_SEED)

production_seed_table <- data.frame(
  replicate = seq_len(N_REP),
  covariate_seed = as.integer(
    MASTER_SEED + 100000L + 1000L * seq_len(N_REP) + 11L
  ),
  sar_seed = as.integer(
    MASTER_SEED + 350000L + 1000L * seq_len(N_REP) + 71L
  ),
  stringsAsFactors = FALSE
)

pilot_seed_table <- data.frame(
  pilot = seq_len(N_PILOT),
  covariate_seed = as.integer(
    MASTER_SEED + 11000000L + 1000L * seq_len(N_PILOT) + 11L
  ),
  sar_seed = as.integer(
    MASTER_SEED + 14000000L + 1000L * seq_len(N_PILOT) + 71L
  ),
  stringsAsFactors = FALSE
)

all_seed_values <- c(
  unlist(production_seed_table[-1L]),
  unlist(pilot_seed_table[-1L]),
  PILOT_MCMC_SEED_BASE + 10000L * seq_len(N_PILOT),
  PRODUCTION_MCMC_SEED_BASE + 10000L * seq_len(N_REP),
  PRODUCTION_PREDICTIVE_SEED_BASE + 10000L * seq_len(N_REP)
)

if (any(all_seed_values > .Machine$integer.max)) {
  stop("At least one deterministic seed exceeds the R integer range.")
}

write.csv(
  production_seed_table,
  file.path(OUT_DIR, "GP_Cov_PU_production_seed_table.csv"),
  row.names = FALSE
)

write.csv(
  pilot_seed_table,
  file.path(OUT_DIR, "GP_Cov_PU_pilot_seed_table.csv"),
  row.names = FALSE
)

## ============================================================================
## 2. INTERNALLY CONSTRUCTED FIXED PRIORS
## ============================================================================

round_prior <- function(x) {
  round(as.numeric(x), PRIOR_ROUND_DIGITS)
}

logit <- function(p) {
  log(p / (1 - p))
}

solve_eta_beta_from_q95 <- function(
    alpha,
    q95_target,
    prob = 0.95
) {
  if (!is.finite(alpha) || alpha <= 0 ||
      !is.finite(q95_target) || q95_target <= 0 || q95_target >= 1) {
    stop("Invalid eta-prior target.")
  }
  
  objective <- function(beta) {
    qbeta(
      prob,
      shape1 = alpha,
      shape2 = beta
    ) - q95_target
  }
  
  lower <- 1e-6
  upper <- 1
  
  while (objective(upper) > 0) {
    upper <- upper * 2
    
    if (upper > 1e6) {
      stop("Could not bracket the eta-prior beta parameter.")
    }
  }
  
  uniroot(
    objective,
    interval = c(lower, upper)
  )$root
}

beta0_low <- logit(BETA0_PREV_RANGE_95[1L])
beta0_high <- logit(BETA0_PREV_RANGE_95[2L])

beta0_prior_for_sampler <- list(
  mean = round_prior(
    0.5 * (beta0_low + beta0_high)
  ),
  sd = round_prior(
    (beta0_high - beta0_low) /
      (2 * qnorm(0.975))
  )
)

eta_beta_unrounded <- solve_eta_beta_from_q95(
  alpha = ETA_ALPHA,
  q95_target = ETA_Q95,
  prob = 0.95
)

eta_prior_for_sampler <- list(
  a_eta = ETA_ALPHA,
  b_eta = round_prior(eta_beta_unrounded)
)

if (!is.finite(beta0_prior_for_sampler$mean) ||
    !is.finite(beta0_prior_for_sampler$sd) ||
    beta0_prior_for_sampler$sd <= 0) {
  stop("The internally constructed beta0 prior is invalid.")
}

if (!is.finite(eta_prior_for_sampler$a_eta) ||
    !is.finite(eta_prior_for_sampler$b_eta) ||
    eta_prior_for_sampler$a_eta <= 0 ||
    eta_prior_for_sampler$b_eta <= 0) {
  stop("The internally constructed eta prior is invalid.")
}

## Every X matrix has ten columns centered and standardized with sample SD = 1.
## Therefore ||X||_F / sqrt(n) is fixed at sqrt(10 * (n - 1) / n).
COVARIATE_MEAN_DIMENSION <- 10L

B_COV_FIXED <- sqrt(
  COVARIATE_MEAN_DIMENSION *
    (N_NODES - 1) /
    N_NODES
)

HALF_T_Q95_UNIT <- qt(
  (0.95 + 1) / 2,
  df = HALFT_DF
)

TAU_SCALE_COV_UNROUNDED <- (
  MEAN_RMS_95_MULT * EFFECTIVE_SD_CENTER_COV
) / (
  HALF_T_Q95_UNIT * B_COV_FIXED
)

TAU_SCALE_COV_ACTIVE <- if (isTRUE(ROUND_ACTIVE_PRIORS)) {
  round_prior(TAU_SCALE_COV_UNROUNDED)
} else {
  TAU_SCALE_COV_UNROUNDED
}

if (!is.finite(TAU_SCALE_COV_ACTIVE) ||
    TAU_SCALE_COV_ACTIVE <= 0) {
  stop("The fixed covariate half-t scale is invalid.")
}

TAU2_COV_DEFAULT <- TAU_SCALE_COV_ACTIVE^2

cat("\n================ INTERNALLY CONSTRUCTED FIXED PRIORS ================\n")

cat(
  "beta0 ~ N(",
  beta0_prior_for_sampler$mean,
  ", ",
  beta0_prior_for_sampler$sd,
  "^2)\n",
  sep = ""
)

cat(
  "eta ~ Beta(",
  eta_prior_for_sampler$a_eta,
  ", ",
  eta_prior_for_sampler$b_eta,
  ")\n",
  sep = ""
)

cat("Single GP effective SD center :", EFFECTIVE_SD_CENTER_COV, "\n")
cat("Effective SD central 95% range: [",
    EFFECTIVE_SD_CENTER_COV / SIGMA_RANGE_FACTOR,
    ", ",
    EFFECTIVE_SD_CENTER_COV * SIGMA_RANGE_FACTOR,
    "]\n",
    sep = "")
cat("Covariate mean dimension      :", COVARIATE_MEAN_DIMENSION, "\n")
cat("Fixed covariate RMS factor    :", B_COV_FIXED, "\n")
cat("tau_cov half-t scale          :", TAU_SCALE_COV_ACTIVE, "\n")
cat("tau_cov initial variance      :", TAU2_COV_DEFAULT, "\n")
cat("Active prior rounding         :", ROUND_ACTIVE_PRIORS, "\n")
cat("Active prior digits           :", PRIOR_ROUND_DIGITS, "\n")

## ============================================================================

## 3. DESIGN AND MCMC CHECKS
## ============================================================================

if (N_PILOT != 5L) {
  stop("This script was requested with exactly five independent pilots.")
}

if (
  N_REP != 100L ||
  N_NODES != 400L ||
  N_TRUE_POS != 45L ||
  N_HIDDEN != 15L ||
  N_OBSERVED_POS != 30L ||
  N_TRUE_ZERO != 355L ||
  N_UNLABELED != 370L
) {
  stop("The fixed latent-truth design has been altered unexpectedly.")
}

if (
  !identical(TRUE_ZERO_IDX, 1:355) ||
  !identical(TRUE_POS_IDX, 356:400)
) {
  stop("The fixed latent-truth row-index design is inconsistent.")
}

if (
  length(COV_ONLY_IDX) != 15L ||
  length(NET_ONLY_IDX) != 15L ||
  length(BOTH_IDX) != 15L ||
  !identical(
    sort(c(COV_ONLY_IDX, NET_ONLY_IDX, BOTH_IDX)),
    TRUE_POS_IDX
  )
) {
  stop("The three mixed-signal positive groups must partition TRUE_POS_IDX.")
}

if (
  length(COVARIATE_SIGNAL_IDX) != 30L ||
  length(NETWORK_SIGNAL_IDX) != 30L ||
  length(intersect(COVARIATE_SIGNAL_IDX, NETWORK_SIGNAL_IDX)) != 15L
) {
  stop("The fixed mixed-signal geometry has been altered unexpectedly.")
}

if (!identical(which(TRUE_T == 1L), TRUE_POS_IDX)) {
  stop("TRUE_T does not match the requested true-positive positions.")
}

if (
  length(LAMBDA_SAR) != 1L ||
  !is.finite(LAMBDA_SAR) ||
  LAMBDA_SAR != 2
) {
  stop("The requested SAR mechanism requires LAMBDA_SAR = 2.")
}

if (
  length(SAR_Z_CLIP) != 1L ||
  !is.finite(SAR_Z_CLIP) ||
  SAR_Z_CLIP <= 0
) {
  stop("SAR_Z_CLIP must be one positive finite number.")
}

if (PILOT_BURNIN <= 0L ||
    PILOT_TAIL_LENGTH <= 0L ||
    PILOT_TAIL_LENGTH > PILOT_BURNIN) {
  stop("Invalid pilot warm-up or tail length.")
}

if (PRODUCTION_BURNIN <= 0L ||
    PRODUCTION_SAMPLING <= 0L ||
    PRODUCTION_THIN <= 0L ||
    PRODUCTION_SAMPLING %% PRODUCTION_THIN != 0L) {
  stop("Invalid production MCMC lengths.")
}

if (!isTRUE(all.equal(PILOT_BLOCKED_FRAC, 1 / 3)) ||
    !isTRUE(all.equal(PRODUCTION_BLOCKED_FRAC, 1 / 3))) {
  stop("Pilot and production blocked fractions must both equal 1/3.")
}

if (length(PREDICTIVE_INTERVAL_LEVELS) != 2L ||
    !isTRUE(all.equal(
      PREDICTIVE_INTERVAL_LEVELS,
      c(0.80, 0.95)
    )) ||
    any(PREDICTIVE_INTERVAL_LEVELS <= 0) ||
    any(PREDICTIVE_INTERVAL_LEVELS >= 1)) {
  stop("PREDICTIVE_INTERVAL_LEVELS must equal c(0.80, 0.95).")
}

if (PREDICTIVE_QUANTILE_TYPE != 1L) {
  stop("Predictive label intervals must use quantile type = 1.")
}

if (length(CORR_PCT_RANGE_ALL) != 2L ||
    any(!is.finite(CORR_PCT_RANGE_ALL)) ||
    CORR_PCT_RANGE_ALL[1L] < 0 ||
    CORR_PCT_RANGE_ALL[2L] > 1 ||
    CORR_PCT_RANGE_ALL[1L] >= CORR_PCT_RANGE_ALL[2L]) {
  stop("CORR_PCT_RANGE_ALL must be an increasing pair inside [0,1].")
}

if (length(ELL_INTERVAL_PROBS) != 2L ||
    length(SIGMA_INTERVAL_PROBS) != 2L ||
    ELL_INTERVAL_PROBS[1L] <= 0 ||
    ELL_INTERVAL_PROBS[2L] >= 1 ||
    SIGMA_INTERVAL_PROBS[1L] <= 0 ||
    SIGMA_INTERVAL_PROBS[2L] >= 1) {
  stop("The ell and sigma interval probabilities are invalid.")
}

if (!identical(as.integer(NU_TAU), as.integer(HALFT_DF))) {
  stop("NU_TAU and HALFT_DF must be identical.")
}

if (isTRUE(FORCE_FRESH_PILOTS)) {
  unlink(PILOT_CHECKPOINT_DIR, recursive = TRUE, force = TRUE)
  dir.create(PILOT_CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(GLOBAL_CALIBRATION_FILE)) {
    file.remove(GLOBAL_CALIBRATION_FILE)
  }
  
  if (file.exists(GLOBAL_CALIBRATION_RDATA)) {
    file.remove(GLOBAL_CALIBRATION_RDATA)
  }
}

if (isTRUE(FORCE_FRESH_PRODUCTION)) {
  unlink(PRODUCTION_CHECKPOINT_DIR, recursive = TRUE, force = TRUE)
  dir.create(
    PRODUCTION_CHECKPOINT_DIR,
    showWarnings = FALSE,
    recursive = TRUE
  )
}

cat("\n================ GP-COV + PU SAR SIMULATION DESIGN ================\n")
cat("Pilot datasets            :", N_PILOT, "\n")
cat("Pilot warm-up             :", PILOT_BURNIN, "\n")
cat("Pilot pooled tail         :", PILOT_TAIL_LENGTH, "per pilot\n")
cat("Production replications   :", N_REP, "\n")
cat("Production burn-in        :", PRODUCTION_BURNIN, "\n")
cat("Production sampling       :", PRODUCTION_SAMPLING, "\n")
cat("Production thinning       :", PRODUCTION_THIN, "\n")
cat("Retained draws/run        :",
    PRODUCTION_SAMPLING / PRODUCTION_THIN, "\n")
cat("Exact GP matrix size      :", N_NODES, "x", N_NODES, "\n")
cat("Predictive interval levels:",
    paste(PREDICTIVE_INTERVAL_LEVELS, collapse = ", "), "\n")
cat("Prior calibration         : exact for every pilot/replication\n")
cat("GP effective SD center    :", EFFECTIVE_SD_CENTER_COV, "\n")
cat("tau half-t df             :", HALFT_DF, "\n")
cat("Application dependency    : none\n")
cat("SAR mechanism             : covariate-driven, lambda =", LAMBDA_SAR, "\n")
cat("SAR z-score clipping      :", SAR_Z_CLIP, "\n")

## ============================================================================
## 4. GENERAL HELPERS
## ============================================================================

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

format_duration <- function(seconds) {
  seconds <- as.numeric(seconds)
  if (!is.finite(seconds)) return(NA_character_)
  hh <- floor(seconds / 3600)
  mm <- floor((seconds - 3600 * hh) / 60)
  ss <- round(seconds - 3600 * hh - 60 * mm)
  sprintf("%02d:%02d:%02d", hh, mm, ss)
}

progress_message <- function(phase, iteration, total) {
  cat(sprintf(
    "[%s] %s %d / %d\n",
    timestamp_now(),
    phase,
    as.integer(iteration),
    as.integer(total)
  ))
}

symmetrize <- function(M) {
  0.5 * (M + t(M))
}

standardize_columns <- function(M) {
  M <- as.matrix(M)
  center <- colMeans(M, na.rm = TRUE)
  scale_value <- apply(M, 2L, stats::sd, na.rm = TRUE)
  center[!is.finite(center)] <- 0
  scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
  scaled <- sweep(M, 2L, center, "-")
  scaled <- sweep(scaled, 2L, scale_value, "/")
  scaled[!is.finite(scaled)] <- 0
  list(scaled = scaled, center = center, scale = scale_value)
}

sqdist <- function(A, B = A) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  aa <- rowSums(A * A)
  bb <- rowSums(B * B)
  pmax(outer(aa, bb, "+") - 2 * tcrossprod(A, B), 0)
}

rbf_kernel_from_D2 <- function(D2, sigma_f = 1, l = 1) {
  sigma_f^2 * exp(-D2 / (2 * l * l))
}

chol_safe <- function(K, jitter = JITTER, max_tries = 8L) {
  K <- symmetrize(as.matrix(K))
  n <- nrow(K)
  for (tt in 0:as.integer(max_tries)) {
    jj <- jitter * (10^tt)
    R <- tryCatch(
      chol(K + diag(jj, n)),
      error = function(e) NULL
    )
    if (!is.null(R)) return(list(R = R, jitter = jj))
  }
  stop("Cholesky factorization failed after jitter escalation.")
}

svd_psd_fallback <- function(K) {
  K <- symmetrize(K)
  E <- tryCatch(eigen(K, symmetric = TRUE), error = function(e) NULL)
  if (!is.null(E)) {
    return(list(U = Re(E$vectors), d = pmax(Re(E$values), 0)))
  }
  out <- svd(K)
  list(U = out$u, d = pmax(out$d, 0))
}

rmvnorm0_chol <- function(Sigma) {
  Sigma <- symmetrize(as.matrix(Sigma))
  p <- nrow(Sigma)
  R <- tryCatch(
    chol(Sigma),
    error = function(e) chol(Sigma + 1e-6 * diag(p))
  )
  as.vector(t(R) %*% rnorm(p))
}

sanitize_cor_2x2 <- function(R, max_abs = 0.95) {
  R <- as.matrix(R)
  if (!all(dim(R) == c(2L, 2L))) return(diag(2))
  rho <- suppressWarnings(as.numeric(R[1L, 2L]))
  if (!is.finite(rho)) rho <- 0
  rho <- max(-abs(max_abs), min(abs(max_abs), rho))
  matrix(c(1, rho, rho, 1), nrow = 2L, byrow = TRUE)
}

estimate_cor_2x2 <- function(
    M,
    fallback = diag(2),
    shrinkage = 0.05,
    max_abs = 0.95,
    min_n = 50L
) {
  M <- as.matrix(M)
  M <- M[complete.cases(M), , drop = FALSE]
  if (nrow(M) < as.integer(min_n) || ncol(M) != 2L) {
    return(sanitize_cor_2x2(fallback, max_abs = max_abs))
  }
  sdv <- apply(M, 2L, sd)
  if (any(!is.finite(sdv)) || any(sdv <= 0)) {
    return(sanitize_cor_2x2(fallback, max_abs = max_abs))
  }
  rho_emp <- suppressWarnings(cor(M[, 1L], M[, 2L]))
  if (!is.finite(rho_emp)) {
    return(sanitize_cor_2x2(fallback, max_abs = max_abs))
  }
  shrinkage <- min(max(as.numeric(shrinkage), 0), 1)
  rho <- (1 - shrinkage) * rho_emp
  sanitize_cor_2x2(
    matrix(c(1, rho, rho, 1), nrow = 2L, byrow = TRUE),
    max_abs = max_abs
  )
}

sanitize_covariance_matrix <- function(S, eig_floor = 1e-10) {
  S <- as.matrix(S)
  if (nrow(S) != ncol(S) || nrow(S) < 1L) {
    stop("S must be a non-empty square matrix.")
  }
  S[!is.finite(S)] <- 0
  ee <- eigen(symmetrize(S), symmetric = TRUE)
  vals <- pmax(Re(ee$values), eig_floor)
  symmetrize(Re(ee$vectors %*% (vals * t(ee$vectors))))
}

safe_cor_from_cov <- function(S) {
  S <- as.matrix(S)
  d <- sqrt(pmax(diag(S), 1e-12))
  C <- S / outer(d, d)
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  C
}

proposal_cov_scaled_2x2 <- function(
    kappa_scalar,
    prior_sd,
    cor_mat,
    jitter = 1e-12
) {
  kappa_scalar <- as.numeric(kappa_scalar)
  prior_sd <- as.numeric(prior_sd)
  if (length(kappa_scalar) != 1L || !is.finite(kappa_scalar) ||
      kappa_scalar <= 0) {
    stop("kappa_scalar must be one positive finite value.")
  }
  if (length(prior_sd) != 2L || any(!is.finite(prior_sd)) ||
      any(prior_sd <= 0)) {
    stop("prior_sd must contain two positive finite values.")
  }
  R <- sanitize_cor_2x2(cor_mat, max_abs = 0.999)
  D <- diag(prior_sd, 2L)
  S <- kappa_scalar * D %*% R %*% D
  symmetrize(S) + diag(jitter, 2L)
}

logN0_eig <- function(x, U, gam, nugget = 0, jitter = 0) {
  nn <- length(x)
  eig <- pmax(as.numeric(gam), 0) + nugget + jitter
  eig <- pmax(eig, .Machine$double.eps)
  ux <- as.vector(crossprod(U, x))
  quad <- sum((ux * ux) / eig)
  ldet <- sum(log(eig))
  -0.5 * (nn * log(2 * pi) + ldet + quad)
}

logN_mu_eig <- function(f, mu, U, gam, nugget = 0, jitter = 0) {
  logN0_eig(f - mu, U, gam, nugget = nugget, jitter = jitter)
}

rInvGamma <- function(n, shape, scale) {
  1 / rgamma(n, shape = shape, rate = scale)
}

with_local_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

validate_probability_vector <- function(probabilities, n_expected, label) {
  probabilities <- as.numeric(probabilities)
  if (length(probabilities) != n_expected) {
    stop(label, " returned the wrong number of probabilities.")
  }
  if (any(!is.finite(probabilities))) {
    stop(label, " returned non-finite probabilities.")
  }
  tolerance <- 1e-8
  if (any(probabilities < -tolerance) ||
      any(probabilities > 1 + tolerance)) {
    stop(label, " returned probabilities outside [0,1].")
  }
  pmin(pmax(probabilities, 0), 1)
}

roc_auc <- function(scores, labels) {
  scores <- as.numeric(scores)
  labels <- as.integer(labels)
  ok <- is.finite(scores) & labels %in% c(0L, 1L)
  scores <- scores[ok]
  labels <- labels[ok]
  n_positive <- sum(labels == 1L)
  n_negative <- sum(labels == 0L)
  if (n_positive == 0L || n_negative == 0L) return(NA_real_)
  ranks_ascending <- rank(scores, ties.method = "average")
  sum_ranks_positive <- sum(ranks_ascending[labels == 1L])
  (
    sum_ranks_positive - n_positive * (n_positive + 1) / 2
  ) / (
    n_positive * n_negative
  )
}

overall_logloss <- function(probabilities, truth, eps = PROBABILITY_EPS) {
  probabilities <- pmin(pmax(as.numeric(probabilities), eps), 1 - eps)
  truth <- as.numeric(truth)
  -mean(
    truth * log(probabilities) +
      (1 - truth) * log(1 - probabilities)
  )
}

hidden_logloss <- function(
    probabilities,
    hidden_idx,
    eps = PROBABILITY_EPS
) {
  probabilities <- pmin(
    pmax(as.numeric(probabilities[hidden_idx]), eps),
    1 - eps
  )
  -mean(log(probabilities))
}

overall_brier <- function(probabilities, truth) {
  mean((as.numeric(probabilities) - as.numeric(truth))^2)
}

hidden_brier <- function(
    probabilities,
    hidden_idx
) {
  mean((as.numeric(probabilities[hidden_idx]) - 1)^2)
}

validate_predictive_label_draws <- function(
    predictive_T_draws,
    n_draws_expected,
    n_nodes_expected,
    method_label
) {
  predictive_T_draws <- as.matrix(predictive_T_draws)
  
  if (!all(dim(predictive_T_draws) == c(
    as.integer(n_draws_expected),
    as.integer(n_nodes_expected)
  ))) {
    stop(
      method_label,
      " returned predictive-label draws with dimensions ",
      paste(dim(predictive_T_draws), collapse = " x "),
      "; expected ",
      n_draws_expected,
      " x ",
      n_nodes_expected,
      "."
    )
  }
  
  if (
    any(!is.finite(predictive_T_draws)) ||
    any(!predictive_T_draws %in% c(0, 1))
  ) {
    stop(method_label, " returned invalid binary predictive-label draws.")
  }
  
  storage.mode(predictive_T_draws) <- "integer"
  predictive_T_draws
}

interval_score_vector <- function(
    L,
    U,
    truth,
    alpha
) {
  L <- as.numeric(L)
  U <- as.numeric(U)
  truth <- as.numeric(truth)
  
  if (
    length(L) != length(U) ||
    length(L) != length(truth)
  ) {
    stop("L, U and truth must have the same length.")
  }
  
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between 0 and 1.")
  }
  
  (U - L) +
    (2 / alpha) * (L - truth) * (truth < L) +
    (2 / alpha) * (truth - U) * (truth > U)
}

binary_predictive_interval <- function(
    predictive_T_draws,
    level,
    quantile_type = PREDICTIVE_QUANTILE_TYPE
) {
  predictive_T_draws <- as.matrix(predictive_T_draws)
  
  if (!is.finite(level) || level <= 0 || level >= 1) {
    stop("level must lie strictly between 0 and 1.")
  }
  
  alpha <- 1 - level
  
  lower <- apply(
    predictive_T_draws,
    2L,
    stats::quantile,
    probs = alpha / 2,
    type = quantile_type,
    names = FALSE
  )
  
  upper <- apply(
    predictive_T_draws,
    2L,
    stats::quantile,
    probs = 1 - alpha / 2,
    type = quantile_type,
    names = FALSE
  )
  
  lower <- as.integer(lower)
  upper <- as.integer(upper)
  
  if (
    any(!lower %in% c(0L, 1L)) ||
    any(!upper %in% c(0L, 1L)) ||
    any(lower > upper)
  ) {
    stop("Invalid empirical interval for binary predictive-label draws.")
  }
  
  list(
    level = level,
    alpha = alpha,
    lower = lower,
    upper = upper
  )
}

predictive_label_metrics <- function(
    predictive_T_draws,
    truth,
    hidden_idx,
    level
) {
  predictive_T_draws <- as.matrix(predictive_T_draws)
  truth <- as.integer(truth)
  hidden_idx <- as.integer(hidden_idx)
  
  interval <- binary_predictive_interval(
    predictive_T_draws = predictive_T_draws,
    level = level
  )
  
  covered <- as.numeric(
    interval$lower <= truth &
      truth <= interval$upper
  )
  
  score <- interval_score_vector(
    L = interval$lower,
    U = interval$upper,
    truth = truth,
    alpha = interval$alpha
  )
  
  list(
    lower = interval$lower,
    upper = interval$upper,
    covered = covered,
    score = score,
    overall_coverage = mean(covered),
    hidden_coverage = mean(covered[hidden_idx]),
    overall_interval_score = mean(score),
    hidden_interval_score = mean(score[hidden_idx])
  )
}

summary_stats <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(c(
      mean = NA_real_,
      sd = NA_real_,
      mc_se = NA_real_,
      q025 = NA_real_,
      median = NA_real_,
      q975 = NA_real_
    ))
  }
  
  sd_x <- if (length(x) > 1L) stats::sd(x) else 0
  
  c(
    mean = mean(x),
    sd = sd_x,
    mc_se = sd_x / sqrt(length(x)),
    q025 = unname(stats::quantile(x, 0.025)),
    median = stats::median(x),
    q975 = unname(stats::quantile(x, 0.975))
  )
}

file_md5_or_na <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path)[1L])
}

median_covariance_2x2 <- function(covariance_list) {
  if (length(covariance_list) < 1L) {
    stop("At least one covariance matrix is required.")
  }
  mats <- lapply(covariance_list, function(S) {
    S <- as.matrix(S)
    if (!all(dim(S) == c(2L, 2L)) || any(!is.finite(S))) {
      stop("All pilot covariance matrices must be finite 2 x 2 matrices.")
    }
    S
  })
  out <- matrix(0, 2L, 2L)
  out[1L, 1L] <- median(vapply(mats, function(S) S[1L, 1L], numeric(1)))
  out[1L, 2L] <- median(vapply(mats, function(S) S[1L, 2L], numeric(1)))
  out[2L, 1L] <- out[1L, 2L]
  out[2L, 2L] <- median(vapply(mats, function(S) S[2L, 2L], numeric(1)))
  sanitize_covariance_matrix(out + COV_JITTER * diag(2L))
}

## ============================================================================

## ============================================================================
## 4A. COVARIATE-DRIVEN SAR LABEL MECHANISM
## ----------------------------------------------------------------------------
## The score below follows the SAR construction and is adapted to
## the mixed-signal DGP. It is an oracle covariate fraud-likeness score:
##
##   larger score  -> covariates look more like the covariate-signal mechanism;
##   smaller score -> covariates look more like the legal/background mechanism.
##
## It is used only to select hidden positives. No fitted method receives the
## score, its block components, its standardized version or its hiding weights.
## ============================================================================

logsumexp2 <- function(a, b) {
  maximum <- pmax(a, b)
  maximum + log(
    exp(a - maximum) +
      exp(b - maximum)
  )
}

log_bivariate_normal_std <- function(
    x,
    y,
    rho
) {
  if (abs(rho) >= 1) {
    stop("rho must lie strictly between -1 and 1.")
  }

  -log(2 * pi) -
    0.5 * log(1 - rho^2) -
    (
      x^2 -
        2 * rho * x * y +
        y^2
    ) / (
      2 * (1 - rho^2)
    )
}

safe_standardize_vector <- function(x) {
  x <- as.numeric(x)
  scale_now <- stats::sd(x)

  if (
    !is.finite(scale_now) ||
    scale_now == 0
  ) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }

  (
    x -
      mean(x)
  ) / scale_now
}

log_trunc_normal_density <- function(
    x,
    mean,
    sd,
    lower,
    upper
) {
  normalizing_probability <-
    stats::pnorm(
      upper,
      mean = mean,
      sd = sd
    ) -
    stats::pnorm(
      lower,
      mean = mean,
      sd = sd
    )

  if (
    !is.finite(normalizing_probability) ||
    normalizing_probability <= 0
  ) {
    stop("Invalid truncated-normal normalizing constant.")
  }

  output <- stats::dnorm(
    x,
    mean = mean,
    sd = sd,
    log = TRUE
  ) -
    log(normalizing_probability)

  output[
    x < lower |
      x > upper
  ] <- -Inf

  output
}

make_oracle_covariate_score <- function(X_raw) {
  X_raw <- as.matrix(X_raw)

  if (
    nrow(X_raw) != N_NODES ||
    ncol(X_raw) != 10L ||
    any(!is.finite(X_raw))
  ) {
    stop(
      "The SAR oracle score requires one finite 400 x 10 raw-covariate matrix."
    )
  }

  ## ------------------------------------------------------------------------
  ## Block 1: exact Gaussian likelihood ratio
  ## ------------------------------------------------------------------------

  score_linear <-
    MU_LINEAR * (
      X_raw[, 1L] +
        X_raw[, 2L]
    ) / SD_LINEAR^2 -
    MU_LINEAR^2 / SD_LINEAR^2

  ## ------------------------------------------------------------------------
  ## Blocks 2 and 4 jointly
  ##
  ## The two covariate-signal profiles share the same latent profile indicator
  ## in the DGP, so their joint positive density is a two-component mixture.
  ## ------------------------------------------------------------------------

  logp_block2_zero <-
    stats::dnorm(
      X_raw[, 3L],
      mean = 0,
      sd = BLOCK2_ZERO_SD_X3,
      log = TRUE
    ) +
    stats::dnorm(
      X_raw[, 4L],
      mean = BLOCK2_X4_CENTER,
      sd = BLOCK2_ZERO_SD_X4,
      log = TRUE
    )

  logp_block2_pos1 <-
    stats::dnorm(
      X_raw[, 3L],
      mean = BLOCK2_X3_OFFSET,
      sd = BLOCK2_POS_SD_X3,
      log = TRUE
    ) +
    stats::dnorm(
      X_raw[, 4L],
      mean = BLOCK2_X4_CENTER -
        BLOCK2_X4_OFFSET,
      sd = BLOCK2_POS_SD_X4,
      log = TRUE
    )

  logp_block2_pos2 <-
    stats::dnorm(
      X_raw[, 3L],
      mean = -BLOCK2_X3_OFFSET,
      sd = BLOCK2_POS_SD_X3,
      log = TRUE
    ) +
    stats::dnorm(
      X_raw[, 4L],
      mean = BLOCK2_X4_CENTER +
        BLOCK2_X4_OFFSET,
      sd = BLOCK2_POS_SD_X4,
      log = TRUE
    )

  parabola_residual <-
    X_raw[, 8L] -
    PARABOLA_COEF *
      X_raw[, 7L]^2

  logp_residual_zero <- stats::dnorm(
    parabola_residual,
    mean = 0,
    sd = SD_PARABOLA_ZERO,
    log = TRUE
  )

  logp_residual_pos1 <- stats::dnorm(
    parabola_residual,
    mean = -DELTA_PARABOLA,
    sd = SD_PARABOLA_POS,
    log = TRUE
  )

  logp_residual_pos2 <- stats::dnorm(
    parabola_residual,
    mean = DELTA_PARABOLA,
    sd = SD_PARABOLA_POS,
    log = TRUE
  )

  logp_block2_parabola_positive <- logsumexp2(
    log(0.5) +
      logp_block2_pos1 +
      logp_residual_pos1,
    log(0.5) +
      logp_block2_pos2 +
      logp_residual_pos2
  )

  logp_block2_parabola_zero <-
    logp_block2_zero +
    logp_residual_zero

  score_profiles_parabola <-
    logp_block2_parabola_positive -
    logp_block2_parabola_zero

  ## ------------------------------------------------------------------------
  ## Block 3: exact disc-versus-ring likelihood ratio
  ## ------------------------------------------------------------------------

  radius <- sqrt(
    X_raw[, 5L]^2 +
      X_raw[, 6L]^2
  )

  logp_disc_zero <- rep(
    -Inf,
    length(radius)
  )

  inside_disc <-
    radius <=
    DISC_RADIUS +
      1e-10

  logp_disc_zero[inside_disc] <-
    -log(
      pi *
        DISC_RADIUS^2
    )

  logp_ring_radius <- log_trunc_normal_density(
    x = radius,
    mean = RING_RADIUS,
    sd = RING_SD,
    lower = 0,
    upper = DISC_RADIUS
  )

  logp_ring_positive <-
    logp_ring_radius -
    log(2 * pi) -
    log(
      pmax(
        radius,
        1e-12
      )
    )

  score_disc_ring <-
    logp_ring_positive -
    logp_disc_zero

  ## ------------------------------------------------------------------------
  ## Block 5: exact opposite-correlation likelihood ratio
  ## ------------------------------------------------------------------------

  logp_corr_zero <- log_bivariate_normal_std(
    X_raw[, 9L],
    X_raw[, 10L],
    rho = RHO_ZERO
  )

  logp_corr_positive <- log_bivariate_normal_std(
    X_raw[, 9L],
    X_raw[, 10L],
    rho = RHO_POS
  )

  score_correlation <-
    logp_corr_positive -
    logp_corr_zero

  block_scores <- cbind(
    linear = score_linear,
    profiles_parabola_joint = score_profiles_parabola,
    disc_ring = score_disc_ring,
    correlation = score_correlation
  )

  if (any(!is.finite(block_scores))) {
    stop("Non-finite values were found in the SAR oracle block scores.")
  }

  oracle_score <- rowSums(
    block_scores
  )

  list(
    score = oracle_score,
    block_scores = block_scores
  )
}

make_replication_signal_type <- function(
    hidden_positive_idx,
    observed_positive_idx
) {
  output <- rep(
    "true zero",
    N_NODES
  )

  output[hidden_positive_idx] <- paste0(
    "hidden: ",
    POSITIVE_SIGNAL_GROUP[hidden_positive_idx]
  )

  output[observed_positive_idx] <- paste0(
    "observed: ",
    POSITIVE_SIGNAL_GROUP[observed_positive_idx]
  )

  output
}

build_sar_label_design <- function(
    X_raw,
    sar_seed
) {
  oracle <- make_oracle_covariate_score(
    X_raw
  )

  positive_scores <- as.numeric(
    oracle$score[
      TRUE_POS_IDX
    ]
  )

  z_cov_positive <- safe_standardize_vector(
    positive_scores
  )

  z_cov_positive <- pmax(
    -SAR_Z_CLIP,
    pmin(
      SAR_Z_CLIP,
      z_cov_positive
    )
  )

  names(z_cov_positive) <- TRUE_POS_IDX

  log_hidden_weights <-
    -LAMBDA_SAR *
    z_cov_positive

  hidden_weights <- exp(
    log_hidden_weights -
      max(log_hidden_weights)
  )

  names(hidden_weights) <- TRUE_POS_IDX

  set.seed(
    as.integer(
      sar_seed
    )
  )

  hidden_positive_idx <- sort(
    sample(
      TRUE_POS_IDX,
      size = N_HIDDEN,
      replace = FALSE,
      prob = hidden_weights
    )
  )

  observed_positive_idx <- sort(
    setdiff(
      TRUE_POS_IDX,
      hidden_positive_idx
    )
  )

  observed_y <- integer(
    N_NODES
  )

  observed_y[
    observed_positive_idx
  ] <- 1L

  signal_type <- make_replication_signal_type(
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = observed_positive_idx
  )

  selection_diagnostics <- data.frame(
    scenario = SCENARIO_NAME,
    node = TRUE_POS_IDX,
    positive_signal_group = POSITIVE_SIGNAL_GROUP[
      TRUE_POS_IDX
    ],
    oracle_covariate_score = positive_scores,
    standardized_oracle_score = as.numeric(
      z_cov_positive
    ),
    hidden_weight = as.numeric(
      hidden_weights
    ),
    hidden_SAR = as.integer(
      TRUE_POS_IDX %in%
        hidden_positive_idx
    ),
    observed_SAR = as.integer(
      TRUE_POS_IDX %in%
        observed_positive_idx
    ),
    stringsAsFactors = FALSE
  )

  selection_summary <- data.frame(
    scenario = SCENARIO_NAME,
    sar_seed = as.integer(sar_seed),
    lambda_SAR = LAMBDA_SAR,
    z_clip = SAR_Z_CLIP,
    n_hidden = length(hidden_positive_idx),
    n_observed_positive = length(observed_positive_idx),
    hidden_covariates_only = sum(
      hidden_positive_idx %in%
        COV_ONLY_IDX
    ),
    hidden_network_only = sum(
      hidden_positive_idx %in%
        NET_ONLY_IDX
    ),
    hidden_both = sum(
      hidden_positive_idx %in%
        BOTH_IDX
    ),
    mean_score_hidden = mean(
      oracle$score[
        hidden_positive_idx
      ]
    ),
    mean_score_observed = mean(
      oracle$score[
        observed_positive_idx
      ]
    ),
    observed_minus_hidden_score = mean(
      oracle$score[
        observed_positive_idx
      ]
    ) -
      mean(
        oracle$score[
          hidden_positive_idx
        ]
      ),
    stringsAsFactors = FALSE
  )

  if (
    length(hidden_positive_idx) != N_HIDDEN ||
    length(observed_positive_idx) != N_OBSERVED_POS ||
    sum(observed_y) != N_OBSERVED_POS ||
    !all(observed_y[hidden_positive_idx] == 0L) ||
    !all(observed_y[observed_positive_idx] == 1L) ||
    !identical(
      sort(c(
        hidden_positive_idx,
        observed_positive_idx
      )),
      TRUE_POS_IDX
    )
  ) {
    stop("The replication-specific SAR label design is inconsistent.")
  }

  list(
    observed_y = observed_y,
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = observed_positive_idx,
    signal_type = signal_type,
    oracle_score = oracle$score,
    oracle_block_scores = oracle$block_scores,
    standardized_positive_score = z_cov_positive,
    hidden_weights = hidden_weights,
    selection_diagnostics = selection_diagnostics,
    selection_summary = selection_summary
  )
}



## ============================================================================
## 5. EXACT GEOMETRY-SPECIFIC COVARIATE-GP PRIOR CALIBRATION
## ============================================================================

fit_lognormal_interval <- function(
    q_low,
    q_high,
    probs = c(0.025, 0.975)
) {
  q_low <- as.numeric(q_low)
  q_high <- as.numeric(q_high)
  
  if (length(q_low) != 1L ||
      length(q_high) != 1L ||
      !is.finite(q_low) ||
      !is.finite(q_high) ||
      q_low <= 0 ||
      q_high <= 0) {
    stop("q_low and q_high must be positive finite scalars.")
  }
  
  qlo <- min(q_low, q_high)
  qhi <- max(q_low, q_high)
  p_low <- probs[1L]
  p_high <- probs[2L]
  
  if (!is.finite(p_low) ||
      !is.finite(p_high) ||
      p_low <= 0 ||
      p_high >= 1 ||
      p_low >= p_high) {
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
  
  pair_d2 <- as.numeric(
    stats::dist(X_std)
  )^2
  
  expected_pairs <- nrow(X_std) *
    (nrow(X_std) - 1) / 2
  
  if (length(pair_d2) != expected_pairs) {
    stop("Unexpected covariate-pair count in prior calibration.")
  }
  
  function(ell) {
    ell <- as.numeric(ell)
    
    if (length(ell) != 1L ||
        !is.finite(ell) ||
        ell <= 0) {
      return(list(
        mean_corr = NA_real_,
        mean_sd = NA_real_
      ))
    }
    
    list(
      mean_corr = mean(
        exp(-pair_d2 / (2 * ell^2))
      ),
      ## The RBF kernel has unit diagonal at unit amplitude.
      mean_sd = 1
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
  
  curve <- curve[
    order(curve$ell),
    ,
    drop = FALSE
  ]
  
  corr_range <- range(
    curve$mean_corr,
    na.rm = TRUE
  )
  
  target_used <- min(
    max(target_corr, corr_range[1L]),
    corr_range[2L]
  )
  
  difference <- curve$mean_corr - target_used
  
  exact_idx <- which(
    abs(difference) <= 1e-12
  )
  
  if (length(exact_idx) > 0L) {
    ell_hat <- curve$ell[exact_idx[1L]]
  } else {
    crossing_idx <- which(
      difference[-nrow(curve)] *
        difference[-1L] <= 0
    )
    
    root_result <- NULL
    
    if (length(crossing_idx) > 0L) {
      crossing_score <- abs(
        difference[crossing_idx]
      ) + abs(
        difference[crossing_idx + 1L]
      )
      
      kk <- crossing_idx[
        which.min(crossing_score)
      ]
      
      log_interval <- log(c(
        curve$ell[kk],
        curve$ell[kk + 1L]
      ))
      
      root_result <- tryCatch(
        uniroot(
          function(log_ell) {
            metric_fun(exp(log_ell))$mean_corr -
              target_used
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
      ordering <- order(
        curve$mean_corr,
        curve$ell
      )
      
      corr_ordered <- curve$mean_corr[ordering]
      ell_ordered <- curve$ell[ordering]
      
      keep <- !duplicated(
        round(corr_ordered, 12)
      )
      
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
        ell_hat <- curve$ell[
          which.min(abs(difference))
        ]
      }
    }
  }
  
  if (!is.finite(ell_hat) ||
      ell_hat <= 0) {
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

calibrate_covariate_gp_prior <- function(
    X_cov,
    calibration_label
) {
  calibration_start <- Sys.time()
  X_cov <- as.matrix(X_cov)
  
  if (nrow(X_cov) != N_NODES ||
      ncol(X_cov) != COVARIATE_MEAN_DIMENSION ||
      any(!is.finite(X_cov))) {
    stop("The covariate prior calibration received an invalid X matrix.")
  }
  
  b_cov_actual <- sqrt(
    sum(X_cov^2) /
      nrow(X_cov)
  )
  
  if (!is.finite(b_cov_actual) ||
      b_cov_actual <= 0) {
    stop("The current covariate RMS factor is invalid.")
  }
  
  if (!isTRUE(all.equal(
    b_cov_actual,
    B_COV_FIXED,
    tolerance = 1e-8
  ))) {
    stop(
      "The covariate design is not standardized as expected: b_cov = ",
      signif(b_cov_actual, 10),
      ", expected ",
      signif(B_COV_FIXED, 10),
      "."
    )
  }
  
  metric_fun <- make_covariate_metric_fun(X_cov)
  
  curve <- build_corr_curve(
    metric_fun = metric_fun,
    ell_bounds = ELL_BOUNDS_COV,
    n_grid = N_ELL_GRID_COARSE
  )
  
  attainable <- range(
    curve$mean_corr,
    na.rm = TRUE
  )
  
  attainable_span <- diff(attainable)
  
  if (!is.finite(attainable_span) ||
      attainable_span <= 1e-10) {
    stop("The attainable covariate correlation range is degenerate.")
  }
  
  targets <- attainable[1L] +
    CORR_PCT_RANGE_ALL * attainable_span
  
  inversion_1 <- invert_corr_curve_refined(
    curve = curve,
    metric_fun = metric_fun,
    target_corr = targets[1L],
    expert_name = "cov"
  )
  
  inversion_2 <- invert_corr_curve_refined(
    curve = curve,
    metric_fun = metric_fun,
    target_corr = targets[2L],
    expert_name = "cov"
  )
  
  ell_low <- min(
    inversion_1$ell,
    inversion_2$ell
  )
  
  ell_high <- max(
    inversion_1$ell,
    inversion_2$ell
  )
  
  prior_ell <- fit_lognormal_interval(
    q_low = ell_low,
    q_high = ell_high,
    probs = ELL_INTERVAL_PROBS
  )
  
  ell_median <- exp(prior_ell$mu)
  unit_kernel_mean_sd <- metric_fun(
    ell_median
  )$mean_sd
  
  if (!is.finite(unit_kernel_mean_sd) ||
      unit_kernel_mean_sd <= 0) {
    stop("Invalid unit-amplitude covariate-GP marginal SD.")
  }
  
  sigma_raw_low <- (
    EFFECTIVE_SD_CENTER_COV /
      SIGMA_RANGE_FACTOR
  ) / unit_kernel_mean_sd
  
  sigma_raw_high <- (
    EFFECTIVE_SD_CENTER_COV *
      SIGMA_RANGE_FACTOR
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
  
  active_hyperprior <- lapply(
    hyperprior_unrounded,
    as.numeric
  )
  
  if (isTRUE(ROUND_ACTIVE_PRIORS)) {
    active_hyperprior <- lapply(
      active_hyperprior,
      round_prior
    )
  }
  
  ## Rounding must never collapse a positive prior SD to zero.
  if (active_hyperprior$log_ell_sd <= 0) {
    active_hyperprior$log_ell_sd <-
      hyperprior_unrounded$log_ell_sd
  }
  
  if (active_hyperprior$log_sig_sd <= 0) {
    active_hyperprior$log_sig_sd <-
      hyperprior_unrounded$log_sig_sd
  }
  
  if (any(!is.finite(unlist(active_hyperprior))) ||
      active_hyperprior$log_ell_sd <= 0 ||
      active_hyperprior$log_sig_sd <= 0) {
    stop("The active geometry-specific covariate prior is invalid.")
  }
  
  theta_prior_mean <- c(
    log_ell_cov = active_hyperprior$log_ell_mean,
    log_sigma_cov = active_hyperprior$log_sig_mean
  )
  
  theta_prior_sd <- c(
    log_ell_cov = active_hyperprior$log_ell_sd,
    log_sigma_cov = active_hyperprior$log_sig_sd
  )
  
  runtime_sec <- as.numeric(difftime(
    Sys.time(),
    calibration_start,
    units = "secs"
  ))
  
  summary_table <- data.frame(
    expert = "cov",
    calibration_label = calibration_label,
    n = nrow(X_cov),
    n_mean_features = ncol(X_cov),
    Z_rms_factor = b_cov_actual,
    Z_rms_factor_from_dimension = B_COV_FIXED,
    
    attainable_corr_min = attainable[1L],
    attainable_corr_max = attainable[2L],
    attainable_corr_span = attainable_span,
    
    target_corr_pct_low = CORR_PCT_RANGE_ALL[1L],
    target_corr_pct_high = CORR_PCT_RANGE_ALL[2L],
    target_corr_low = targets[1L],
    target_corr_high = targets[2L],
    
    achieved_corr_at_first_target =
      inversion_1$achieved_corr,
    achieved_corr_at_second_target =
      inversion_2$achieved_corr,
    
    ell_interval_low = ell_low,
    ell_interval_high = ell_high,
    log_ell_mean_unrounded = prior_ell$mu,
    log_ell_sd_unrounded = prior_ell$tau,
    ell_median = prior_ell$median,
    
    unit_kernel_mean_marginal_sd =
      unit_kernel_mean_sd,
    
    effective_sd_target_low =
      EFFECTIVE_SD_CENTER_COV / SIGMA_RANGE_FACTOR,
    effective_sd_target_center =
      EFFECTIVE_SD_CENTER_COV,
    effective_sd_target_high =
      EFFECTIVE_SD_CENTER_COV * SIGMA_RANGE_FACTOR,
    
    sigma_raw_interval_low = sigma_raw_low,
    sigma_raw_interval_high = sigma_raw_high,
    log_sig_mean_unrounded = prior_sigma$mu,
    log_sig_sd_unrounded = prior_sigma$tau,
    sigma_raw_median = prior_sigma$median,
    effective_sigma_median_check =
      prior_sigma$median * unit_kernel_mean_sd,
    
    log_ell_mean_active =
      active_hyperprior$log_ell_mean,
    log_ell_sd_active =
      active_hyperprior$log_ell_sd,
    log_sig_mean_active =
      active_hyperprior$log_sig_mean,
    log_sig_sd_active =
      active_hyperprior$log_sig_sd,
    
    half_t_df = HALFT_DF,
    mean_rms_95_mult = MEAN_RMS_95_MULT,
    tau_scale_unrounded =
      TAU_SCALE_COV_UNROUNDED,
    tau_scale_active =
      TAU_SCALE_COV_ACTIVE,
    tau2_default =
      TAU2_COV_DEFAULT,
    
    active_prior_rounding =
      ROUND_ACTIVE_PRIORS,
    active_prior_digits =
      PRIOR_ROUND_DIGITS,
    prior_calibration_runtime_sec =
      runtime_sec,
    
    stringsAsFactors = FALSE
  )
  
  list(
    calibration_label = calibration_label,
    active_hyperprior = active_hyperprior,
    hyperprior_unrounded = hyperprior_unrounded,
    tau_scale = TAU_SCALE_COV_ACTIVE,
    tau_scale_unrounded = TAU_SCALE_COV_UNROUNDED,
    tau2_default = TAU2_COV_DEFAULT,
    theta_prior_mean = theta_prior_mean,
    theta_prior_sd = theta_prior_sd,
    summary = summary_table,
    curve = curve,
    runtime_sec = runtime_sec
  )
}

## ============================================================================
## 6. COVARIATE DATA-GENERATING PROCESS
## ============================================================================

rmvnorm2 <- function(
    n,
    rho,
    sd1 = 1,
    sd2 = 1,
    mean = c(0, 0)
) {
  if (abs(rho) >= 1) {
    stop(
      "rho must lie strictly between -1 and 1."
    )
  }
  
  Sigma <- matrix(
    c(
      sd1^2,
      rho * sd1 * sd2,
      rho * sd1 * sd2,
      sd2^2
    ),
    nrow = 2L,
    byrow = TRUE
  )
  
  Z <- matrix(
    rnorm(
      2L * n
    ),
    ncol = 2L
  )
  
  X <- Z %*% chol(
    Sigma
  )
  
  X[, 1L] <- X[, 1L] +
    mean[1L]
  
  X[, 2L] <- X[, 2L] +
    mean[2L]
  
  X
}

rtrunc_normal <- function(
    n,
    mean,
    sd,
    lower,
    upper
) {
  n <- as.integer(n)
  
  if (n <= 0L) {
    return(
      numeric(0)
    )
  }
  
  if (
    !is.finite(mean) ||
    !is.finite(sd) ||
    sd <= 0 ||
    !is.finite(lower) ||
    !is.finite(upper) ||
    lower >= upper
  ) {
    stop(
      "Invalid truncated-normal parameters."
    )
  }
  
  out <- numeric(0)
  
  while (length(out) < n) {
    remaining <- n -
      length(out)
    
    candidate <- rnorm(
      max(
        4L * remaining,
        50L
      ),
      mean = mean,
      sd = sd
    )
    
    candidate <- candidate[
      candidate >= lower &
        candidate <= upper
    ]
    
    out <- c(
      out,
      candidate
    )
  }
  
  out[
    seq_len(n)
  ]
}

make_covariate_set <- function(
    seed,
    signal_idx = COVARIATE_SIGNAL_IDX
) {
  set.seed(seed)
  
  signal_idx <- sort(
    unique(
      as.integer(signal_idx)
    )
  )
  
  if (length(signal_idx) != 30L) {
    stop(
      "The covariate DGP requires exactly 30 signal nodes."
    )
  }
  
  if (
    any(signal_idx < 1L) ||
    any(signal_idx > N_NODES)
  ) {
    stop(
      "signal_idx contains invalid node positions."
    )
  }
  
  background_idx <- setdiff(
    seq_len(N_NODES),
    signal_idx
  )
  
  ## Randomly allocate the 30 signal nodes to the two positive profiles,
  ## retaining exactly 15 nodes in each. This avoids confounding a profile
  ## deterministically with observed/hidden status.
  profile_order <- sample(
    signal_idx,
    length(signal_idx),
    replace = FALSE
  )
  
  profile_1_idx <- sort(
    profile_order[1:15]
  )
  
  profile_2_idx <- sort(
    profile_order[16:30]
  )
  
  X <- matrix(
    NA_real_,
    nrow = N_NODES,
    ncol = 10L,
    dimnames = list(
      NULL,
      paste0(
        "X",
        seq_len(10L)
      )
    )
  )
  
  n_background <- length(
    background_idx
  )
  
  n_signal <- length(
    signal_idx
  )
  
  n_profile_1 <- length(
    profile_1_idx
  )
  
  n_profile_2 <- length(
    profile_2_idx
  )
  
  ## Block 1: mild linear signal.
  X[background_idx, 1L] <- rnorm(
    n_background,
    mean = 0,
    sd = SD_LINEAR
  )
  
  X[background_idx, 2L] <- rnorm(
    n_background,
    mean = 0,
    sd = SD_LINEAR
  )
  
  X[signal_idx, 1L] <- rnorm(
    n_signal,
    mean = MU_LINEAR,
    sd = SD_LINEAR
  )
  
  X[signal_idx, 2L] <- rnorm(
    n_signal,
    mean = MU_LINEAR,
    sd = SD_LINEAR
  )
  
  ## Block 2: broad background cloud and two signal profiles.
  X[background_idx, 3L] <- rnorm(
    n_background,
    mean = 0,
    sd = BLOCK2_ZERO_SD_X3
  )
  
  X[background_idx, 4L] <- rnorm(
    n_background,
    mean = BLOCK2_X4_CENTER,
    sd = BLOCK2_ZERO_SD_X4
  )
  
  X[profile_1_idx, 3L] <- rnorm(
    n_profile_1,
    mean = BLOCK2_X3_OFFSET,
    sd = BLOCK2_POS_SD_X3
  )
  
  X[profile_1_idx, 4L] <- rnorm(
    n_profile_1,
    mean = BLOCK2_X4_CENTER -
      BLOCK2_X4_OFFSET,
    sd = BLOCK2_POS_SD_X4
  )
  
  X[profile_2_idx, 3L] <- rnorm(
    n_profile_2,
    mean = -BLOCK2_X3_OFFSET,
    sd = BLOCK2_POS_SD_X3
  )
  
  X[profile_2_idx, 4L] <- rnorm(
    n_profile_2,
    mean = BLOCK2_X4_CENTER +
      BLOCK2_X4_OFFSET,
    sd = BLOCK2_POS_SD_X4
  )
  
  ## Block 3: background disc and signal ring.
  theta_background <- runif(
    n_background,
    min = 0,
    max = 2 * pi
  )
  
  radius_background <- DISC_RADIUS *
    sqrt(
      runif(
        n_background
      )
    )
  
  X[background_idx, 5L] <- radius_background *
    cos(
      theta_background
    )
  
  X[background_idx, 6L] <- radius_background *
    sin(
      theta_background
    )
  
  theta_signal <- runif(
    n_signal,
    min = 0,
    max = 2 * pi
  )
  
  radius_signal <- rtrunc_normal(
    n = n_signal,
    mean = RING_RADIUS,
    sd = RING_SD,
    lower = 0,
    upper = DISC_RADIUS
  )
  
  X[signal_idx, 5L] <- radius_signal *
    cos(
      theta_signal
    )
  
  X[signal_idx, 6L] <- radius_signal *
    sin(
      theta_signal
    )
  
  ## Block 4: parabolic residual structure.
  U <- rnorm(
    N_NODES
  )
  
  X[, 7L] <- U
  
  expected_x8 <- PARABOLA_COEF *
    U^2
  
  X[background_idx, 8L] <- expected_x8[
    background_idx
  ] +
    rnorm(
      n_background,
      mean = 0,
      sd = SD_PARABOLA_ZERO
    )
  
  X[profile_1_idx, 8L] <- expected_x8[
    profile_1_idx
  ] -
    DELTA_PARABOLA +
    rnorm(
      n_profile_1,
      mean = 0,
      sd = SD_PARABOLA_POS
    )
  
  X[profile_2_idx, 8L] <- expected_x8[
    profile_2_idx
  ] +
    DELTA_PARABOLA +
    rnorm(
      n_profile_2,
      mean = 0,
      sd = SD_PARABOLA_POS
    )
  
  ## Block 5: opposite correlation structure.
  X[background_idx, 9L:10L] <- rmvnorm2(
    n = n_background,
    rho = RHO_ZERO
  )
  
  X[signal_idx, 9L:10L] <- rmvnorm2(
    n = n_signal,
    rho = RHO_POS
  )
  
  standardized <- standardize_columns(
    X
  )
  
  covariate_profile <- rep(
    "background",
    N_NODES
  )
  
  covariate_profile[
    profile_1_idx
  ] <- "signal profile 1"
  
  covariate_profile[
    profile_2_idx
  ] <- "signal profile 2"
  
  list(
    seed = seed,
    X_raw = X,
    X_std = standardized$scaled,
    center = standardized$center,
    scale = standardized$scale,
    signal_idx = signal_idx,
    background_idx = background_idx,
    profile_1_idx = profile_1_idx,
    profile_2_idx = profile_2_idx,
    covariate_profile = covariate_profile
  )
}

## ============================================================================

## ============================================================================
## 7. PU LIKELIHOOD AND ETA SLICE SAMPLER
## ============================================================================

pu_loglik_vec <- function(f, y, eta, eps = 1e-12) {
  p_latent <- plogis(as.numeric(f))
  p_observed <- (1 - eta) * p_latent
  p_observed <- pmin(pmax(p_observed, eps), 1 - eps)
  y <- as.numeric(y)
  
  sum(
    y * log(p_observed) +
      (1 - y) * log(1 - p_observed)
  )
}

pu_gradloglik_vec <- function(f, y, eta, eps = 1e-12) {
  p_latent <- plogis(as.numeric(f))
  p_observed <- (1 - eta) * p_latent
  p_observed <- pmin(pmax(p_observed, eps), 1 - eps)
  y <- as.numeric(y)
  
  (1 - eta) * p_latent * (1 - p_latent) *
    (
      y / p_observed -
        (1 - y) / (1 - p_observed)
    )
}

eta_logpost_u <- function(
    u,
    f_fixed,
    Y,
    a_eta,
    b_eta,
    eps = 1e-12
) {
  eta <- plogis(u)
  eta <- pmin(pmax(eta, eps), 1 - eps)
  
  pu_loglik_vec(
    f_fixed,
    Y,
    eta,
    eps = eps
  ) +
    dbeta(
      eta,
      a_eta,
      b_eta,
      log = TRUE
    ) +
    log(eta) +
    log(1 - eta)
}

slice_eta_step <- function(
    eta,
    f_fixed,
    Y,
    a_eta,
    b_eta,
    width = 0.50,
    m = 100L,
    eps = 1e-12
) {
  u0 <- qlogis(
    pmin(pmax(eta, eps), 1 - eps)
  )
  
  log_density_0 <- eta_logpost_u(
    u0,
    f_fixed,
    Y,
    a_eta,
    b_eta,
    eps = eps
  )
  
  log_height <- log_density_0 - rexp(1)
  
  offset <- runif(1, 0, width)
  left <- u0 - offset
  right <- left + width
  
  J <- floor(runif(1, 0, m))
  K <- (m - 1L) - J
  n_evaluations <- 1L
  
  while (J > 0L) {
    value_left <- eta_logpost_u(
      left,
      f_fixed,
      Y,
      a_eta,
      b_eta,
      eps = eps
    )
    
    n_evaluations <- n_evaluations + 1L
    
    if (!is.finite(value_left) ||
        value_left <= log_height) {
      break
    }
    
    left <- left - width
    J <- J - 1L
  }
  
  while (K > 0L) {
    value_right <- eta_logpost_u(
      right,
      f_fixed,
      Y,
      a_eta,
      b_eta,
      eps = eps
    )
    
    n_evaluations <- n_evaluations + 1L
    
    if (!is.finite(value_right) ||
        value_right <= log_height) {
      break
    }
    
    right <- right + width
    K <- K - 1L
  }
  
  repeat {
    u1 <- runif(1, left, right)
    
    value_1 <- eta_logpost_u(
      u1,
      f_fixed,
      Y,
      a_eta,
      b_eta,
      eps = eps
    )
    
    n_evaluations <- n_evaluations + 1L
    
    if (is.finite(value_1) &&
        value_1 >= log_height) {
      break
    }
    
    if (u1 < u0) {
      left <- u1
    } else {
      right <- u1
    }
  }
  
  eta_new <- plogis(u1)
  eta_new <- pmin(pmax(eta_new, eps), 1 - eps)
  
  list(
    eta = eta_new,
    u = u1,
    moved = as.numeric(abs(u1 - u0) > 1e-12),
    n_eval = n_evaluations
  )
}

hidden_probability_given_zero <- function(
    p_latent,
    eta
) {
  p_latent <- pmin(
    pmax(as.numeric(p_latent), 0),
    1
  )
  
  eta <- pmin(
    pmax(as.numeric(eta), ETA_EPS),
    1 - ETA_EPS
  )
  
  denominator <- 1 - (1 - eta) * p_latent
  denominator <- pmax(
    denominator,
    .Machine$double.eps
  )
  
  output <- eta * p_latent / denominator
  
  pmin(pmax(output, 0), 1)
}

## ============================================================================

## 8. EXACT GP-COV + PU CHAIN
## ============================================================================

run_gp_cov_pu_exact_chain <- function(
    Y,
    X_cov,
    hyperprior_cov,
    m_beta,
    s_beta,
    tau2_default,
    scale_tau,
    a_eta,
    b_eta,
    n_burnin,
    n_sample = 0L,
    thin = 1L,
    blocked_frac = 1 / 3,
    min_blocked_warmup = 5000L,
    init_state = NULL,
    delta_init = 3.5,
    kappa_cov_init = 1,
    proposal_cor_init = diag(2),
    adapt_block_cor = TRUE,
    cor_adapt_start = 2000L,
    cor_adapt_interval = 250L,
    cor_adapt_window = 5000L,
    cor_shrinkage = 0.05,
    cor_max_abs = 0.95,
    joint_cov_base_init = NULL,
    joint_theta_scale_init = 1,
    f_refresh_after_joint = TRUE,
    eta_slice_width = 0.50,
    initial_eta = 0.20,
    nu_tau = 3,
    sample_tau = TRUE,
    adapt = TRUE,
    adapt_start = 2000L,
    adapt_interval = 25L,
    target_f = 0.55,
    target_theta = 0.25,
    target_joint = 0.25,
    adapt_rate_delta = 0.05,
    adapt_rate_theta = 0.05,
    adapt_rate_joint_theta = 0.5,
    t0_adapt = 10,
    power_adapt = 0.60,
    delta_min = 1e-4,
    delta_max = 50,
    kappa_min = 1e-8,
    kappa_max = 100,
    joint_theta_scale_min = 1e-4,
    joint_theta_scale_max = 20,
    emp_cov_tail = 5000L,
    cov_jitter = 1e-6,
    max_log_sigma_cov = log(5),
    min_log_sigma = log(1e-3),
    logl_cov_bounds = c(log(1e-3), log(1e3)),
    jitter = 1e-8,
    state_tail_length = 0L,
    predictive_seed = NULL,
    phase_label = "GP-Cov + PU",
    burnin_progress_every = 5000L,
    sample_progress_every = 5000L,
    verbose = TRUE
) {
  Y <- as.numeric(Y)
  X_cov <- as.matrix(X_cov)
  n <- length(Y)
  p_mean <- ncol(X_cov)
  
  n_burnin <- as.integer(n_burnin)
  n_sample <- as.integer(n_sample)
  thin <- as.integer(thin)
  min_blocked_warmup <- as.integer(min_blocked_warmup)
  adapt_start <- as.integer(adapt_start)
  adapt_interval <- as.integer(adapt_interval)
  cor_adapt_start <- as.integer(cor_adapt_start)
  cor_adapt_interval <- as.integer(cor_adapt_interval)
  cor_adapt_window <- as.integer(cor_adapt_window)
  emp_cov_tail <- as.integer(emp_cov_tail)
  state_tail_length <- as.integer(state_tail_length)
  
  if (!all(Y %in% c(0, 1))) stop("Y must contain only 0 and 1.")
  if (nrow(X_cov) != n) stop("X_cov and Y have incompatible dimensions.")
  if (p_mean != 10L) {
    stop("The GP-Cov + PU simulation expects exactly ten covariate columns.")
  }
  if (n_burnin < 0L || n_sample < 0L || thin <= 0L) {
    stop("Invalid chain lengths.")
  }
  if (n_sample > 0L && n_sample %% thin != 0L) {
    stop("n_sample must be divisible by thin.")
  }
  if (state_tail_length < 0L || state_tail_length > n_burnin) {
    stop("state_tail_length must lie between zero and n_burnin.")
  }
  if (!is.finite(scale_tau) || scale_tau <= 0 ||
      !is.finite(tau2_default) || tau2_default <= 0) {
    stop("The tau prior settings must be positive and finite.")
  }
  if (!is.finite(a_eta) || a_eta <= 0 ||
      !is.finite(b_eta) || b_eta <= 0) {
    stop("The eta prior settings must be positive and finite.")
  }
  
  D2 <- sqdist(X_cov)
  hp <- hyperprior_cov
  prior_sd_theta <- pmax(c(hp$log_ell_sd, hp$log_sig_sd), 1e-8)
  names(prior_sd_theta) <- c("log_ell_cov", "log_sigma_cov")
  
  blocked_nominal <- as.integer(floor(n_burnin * blocked_frac))
  blocked_iters <- if (n_burnin > 0L) {
    min(n_burnin, max(blocked_nominal, min_blocked_warmup))
  } else {
    0L
  }
  joint_tune_iters <- n_burnin - blocked_iters
  n_saved_expected <- if (n_sample > 0L) n_sample %/% thin else 0L
  
  normalize_scalar <- function(x, arg_name) {
    x <- as.numeric(x)
    if (length(x) != 1L || !is.finite(x) || x <= 0) {
      stop(arg_name, " must be one positive finite scalar.")
    }
    x
  }
  
  theta_ok <- function(theta) {
    all(is.finite(theta)) &&
      theta[1L] >= logl_cov_bounds[1L] &&
      theta[1L] <= logl_cov_bounds[2L] &&
      theta[2L] >= min_log_sigma &&
      theta[2L] <= max_log_sigma_cov
  }
  
  logprior_theta <- function(theta) {
    if (!theta_ok(theta)) return(-Inf)
    dnorm(
      theta[1L],
      mean = hp$log_ell_mean,
      sd = max(hp$log_ell_sd, 1e-8),
      log = TRUE
    ) +
      dnorm(
        theta[2L],
        mean = hp$log_sig_mean,
        sd = max(hp$log_sig_sd, 1e-8),
        log = TRUE
      )
  }
  
  K_cov_from_theta <- function(theta) {
    rbf_kernel_from_D2(
      D2,
      sigma_f = exp(theta[2L]),
      l = exp(theta[1L])
    )
  }
  
  Csum_svd <- function(theta) {
    K <- symmetrize(K_cov_from_theta(theta)) + diag(jitter, n)
    svd_psd_fallback(K)
  }
  
  Z_mean <- cbind(beta0 = 1, X_cov)
  feature_names <- colnames(X_cov)
  if (is.null(feature_names)) feature_names <- paste0("X", seq_len(p_mean))
  colnames(Z_mean) <- c("beta0", paste0("gamma_", feature_names))
  p_coef <- ncol(Z_mean)
  prior_mean_coef <- c(m_beta, rep(0, p_mean))
  names(prior_mean_coef) <- colnames(Z_mean)
  
  if (is.null(init_state)) {
    theta_cov <- c(hp$log_ell_mean, hp$log_sig_mean)
    coef_vec <- prior_mean_coef
    tau2_cov <- as.numeric(tau2_default)
    aux_tau_cov <- rInvGamma(
      1L,
      shape = 1 / 2,
      scale = 1 / (scale_tau^2)
    )
    eta <- pmin(
      pmax(initial_eta, ETA_EPS),
      1 - ETA_EPS
    )
    f_init_supplied <- NULL
  } else {
    theta_cov <- as.numeric(init_state$theta_log)
    coef_vec <- as.numeric(init_state$beta)
    tau2_cov <- as.numeric(init_state$tau2_cov)
    aux_tau_cov <- as.numeric(init_state$aux_tau_cov)
    eta <- as.numeric(init_state$eta)
    
    if (length(eta) != 1L ||
        !is.finite(eta)) {
      stop("The initial eta state is invalid.")
    }
    
    eta <- pmin(
      pmax(eta, ETA_EPS),
      1 - ETA_EPS
    )
    
    f_init_supplied <- init_state$f
  }
  
  if (length(theta_cov) != 2L || any(!is.finite(theta_cov))) {
    stop("The initial theta state is invalid.")
  }
  theta_cov[1L] <- min(
    max(theta_cov[1L], logl_cov_bounds[1L]),
    logl_cov_bounds[2L]
  )
  theta_cov[2L] <- min(
    max(theta_cov[2L], min_log_sigma),
    max_log_sigma_cov
  )
  if (!theta_ok(theta_cov)) stop("The initial theta state violates bounds.")
  
  if (length(coef_vec) != p_coef || any(!is.finite(coef_vec))) {
    stop("The initial beta/gamma state has the wrong length or invalid values.")
  }
  names(coef_vec) <- colnames(Z_mean)
  if (!is.finite(tau2_cov) || tau2_cov <= 0 ||
      !is.finite(aux_tau_cov) || aux_tau_cov <= 0) {
    stop("The initial tau hierarchy state is invalid.")
  }
  
  mu_vec <- as.vector(Z_mean %*% coef_vec)
  if (is.null(f_init_supplied)) {
    f <- mu_vec
  } else {
    f <- as.numeric(f_init_supplied)
    if (length(f) != n || any(!is.finite(f))) {
      stop("The supplied initial f state is invalid.")
    }
  }
  
  svd_sum <- Csum_svd(theta_cov)
  U_sum <- svd_sum$U
  gam_sum <- pmax(svd_sum$d, 0)
  
  lf <- pu_loglik_vec(f, Y, eta)
  gf <- pu_gradloglik_vec(f, Y, eta)
  
  kappa_cov <- normalize_scalar(kappa_cov_init, "kappa_cov_init")
  proposal_cor <- sanitize_cor_2x2(
    proposal_cor_init,
    max_abs = cor_max_abs
  )
  joint_theta_scale <- normalize_scalar(
    joint_theta_scale_init,
    "joint_theta_scale_init"
  )
  delta_state <- normalize_scalar(delta_init, "delta_init")
  
  sample_beta_gamma_gibbs <- function(f, U, gam, tau2_cov) {
    prior_var <- c(s_beta^2, rep(tau2_cov, p_mean))
    prior_var <- pmax(prior_var, .Machine$double.eps)
    eig <- pmax(as.numeric(gam), .Machine$double.eps)
    UtZ <- crossprod(U, Z_mean)
    KinvZ <- U %*% sweep(UtZ, 1L, eig, "/")
    Utf <- as.vector(crossprod(U, f))
    Kinvf <- as.vector(U %*% (Utf / eig))
    Q <- symmetrize(
      crossprod(Z_mean, KinvZ) + diag(1 / prior_var, p_coef)
    )
    h <- as.vector(crossprod(Z_mean, Kinvf)) +
      prior_mean_coef / prior_var
    R <- chol_safe(Q, jitter = 1e-10, max_tries = 8L)$R
    mean_coef <- as.vector(backsolve(R, forwardsolve(t(R), h)))
    mean_coef + as.vector(backsolve(R, rnorm(p_coef)))
  }
  
  update_coef_tau <- function(
    f,
    U,
    gam,
    coef_vec,
    tau2_cov,
    aux_tau_cov
  ) {
    coef_new <- sample_beta_gamma_gibbs(f, U, gam, tau2_cov)
    names(coef_new) <- colnames(Z_mean)
    if (isTRUE(sample_tau)) {
      gamma_new <- coef_new[-1L]
      tau2_cov <- rInvGamma(
        1L,
        shape = (nu_tau + p_mean) / 2,
        scale = (nu_tau / aux_tau_cov) + sum(gamma_new^2) / 2
      )
      aux_tau_cov <- rInvGamma(
        1L,
        shape = (nu_tau + 1) / 2,
        scale = (1 / scale_tau^2) + (nu_tau / tau2_cov)
      )
    }
    list(
      coef_vec = coef_new,
      tau2_cov = tau2_cov,
      aux_tau_cov = aux_tau_cov,
      mu_vec = as.vector(Z_mean %*% coef_new)
    )
  }
  
  update_eta_coef_tau <- function(
    f,
    eta_current,
    U,
    gam,
    coef_vec,
    tau2_cov,
    aux_tau_cov
  ) {
    eta_update <- slice_eta_step(
      eta = eta_current,
      f_fixed = f,
      Y = Y,
      a_eta = a_eta,
      b_eta = b_eta,
      width = eta_slice_width,
      m = 100L
    )
    
    eta_new <- eta_update$eta
    
    coef_update <- update_coef_tau(
      f = f,
      U = U,
      gam = gam,
      coef_vec = coef_vec,
      tau2_cov = tau2_cov,
      aux_tau_cov = aux_tau_cov
    )
    
    list(
      eta = eta_new,
      lf = pu_loglik_vec(f, Y, eta_new),
      gf = pu_gradloglik_vec(f, Y, eta_new),
      coef_vec = coef_update$coef_vec,
      tau2_cov = coef_update$tau2_cov,
      aux_tau_cov = coef_update$aux_tau_cov,
      mu_vec = coef_update$mu_vec,
      eta_moved = eta_update$moved,
      eta_n_eval = eta_update$n_eval
    )
  }
  
  adapt_log_step <- function(
    step,
    acc,
    target,
    iter,
    rate,
    t0,
    power,
    min_step,
    max_step
  ) {
    gain <- rate / ((iter + t0)^power)
    step_new <- step * exp(gain * (acc - target))
    max(min_step, min(max_step, step_new))
  }
  
  agrad_f_step <- function(f, mu, lf, gf, eta, U, gam, delta) {
    nn <- length(f)
    g <- f - mu
    z <- g + (delta / 2) * gf + sqrt(delta / 2) * rnorm(nn)
    sqrt_lt <- sqrt(
      (pmax(gam, 0) * delta) /
        pmax(delta + 2 * pmax(gam, 0), .Machine$double.eps)
    )
    sqrt_lt[!is.finite(sqrt_lt)] <- 0
    temp1 <- crossprod(U, (2 / delta) * z)
    temp2 <- sqrt_lt * temp1 + rnorm(nn)
    temp3 <- sqrt_lt * temp2
    g_prop <- as.vector(U %*% temp3)
    f_prop <- mu + g_prop
    lf_prop <- pu_loglik_vec(f_prop, Y, eta)
    gf_prop <- pu_gradloglik_vec(f_prop, Y, eta)
    gzy <- as.numeric(crossprod(
      z - g_prop - (delta / 4) * gf_prop,
      gf_prop
    ))
    gzx <- as.numeric(crossprod(
      z - g - (delta / 4) * gf,
      gf
    ))
    log_alpha <- (lf_prop - lf) + (gzy - gzx)
    acc <- log(runif(1)) < min(0, log_alpha)
    if (acc) {
      list(f = f_prop, lf = lf_prop, gf = gf_prop, acc = TRUE)
    } else {
      list(f = f, lf = lf, gf = gf, acc = FALSE)
    }
  }
  
  theta_block_step <- function(
    theta,
    f_fixed,
    mu_fixed,
    U_curr,
    gam_curr,
    kappa_scalar,
    cor_mat
  ) {
    lp_curr <- logprior_theta(theta)
    logN_curr <- logN_mu_eig(f_fixed, mu_fixed, U_curr, gam_curr)
    Sigma_block <- proposal_cov_scaled_2x2(
      kappa_scalar = kappa_scalar,
      prior_sd = prior_sd_theta,
      cor_mat = cor_mat
    )
    theta_prop <- theta + rmvnorm0_chol(Sigma_block)
    lp_prop <- logprior_theta(theta_prop)
    if (!is.finite(lp_prop)) {
      return(list(theta = theta, U = U_curr, gam = gam_curr, acc = FALSE))
    }
    svd_prop <- Csum_svd(theta_prop)
    U_prop <- svd_prop$U
    gam_prop <- pmax(svd_prop$d, 0)
    logN_prop <- logN_mu_eig(f_fixed, mu_fixed, U_prop, gam_prop)
    log_alpha <- (logN_prop - logN_curr) + (lp_prop - lp_curr)
    acc <- log(runif(1)) < min(0, log_alpha)
    if (acc) {
      list(theta = theta_prop, U = U_prop, gam = gam_prop, acc = TRUE)
    } else {
      list(theta = theta, U = U_curr, gam = gam_curr, acc = FALSE)
    }
  }
  
  joint_step <- function(
    f,
    mu,
    lf,
    gf,
    theta,
    eta,
    U_curr,
    gam_curr,
    delta,
    Sigma_theta
  ) {
    nn <- length(f)
    g <- f - mu
    z <- g + (delta / 2) * gf + sqrt(delta / 2) * rnorm(nn)
    theta_prop <- theta + rmvnorm0_chol(Sigma_theta)
    lp_curr <- logprior_theta(theta)
    lp_prop <- logprior_theta(theta_prop)
    if (!is.finite(lp_prop)) {
      return(list(
        f = f,
        lf = lf,
        gf = gf,
        theta = theta,
        U = U_curr,
        gam = gam_curr,
        acc = FALSE
      ))
    }
    svd_prop <- Csum_svd(theta_prop)
    U_prop <- svd_prop$U
    gam_prop <- pmax(svd_prop$d, 0)
    sqrt_lt_prop <- sqrt(
      (pmax(gam_prop, 0) * delta) /
        pmax(delta + 2 * pmax(gam_prop, 0), .Machine$double.eps)
    )
    sqrt_lt_prop[!is.finite(sqrt_lt_prop)] <- 0
    temp1 <- crossprod(U_prop, (2 / delta) * z)
    temp2 <- sqrt_lt_prop * temp1 + rnorm(nn)
    temp3 <- sqrt_lt_prop * temp2
    g_prop <- as.vector(U_prop %*% temp3)
    f_prop <- mu + g_prop
    lf_prop <- pu_loglik_vec(f_prop, Y, eta)
    gf_prop <- pu_gradloglik_vec(f_prop, Y, eta)
    gzy <- as.numeric(crossprod(
      z - g_prop - (delta / 4) * gf_prop,
      gf_prop
    ))
    gzx <- as.numeric(crossprod(
      z - g - (delta / 4) * gf,
      gf
    ))
    logZ_curr <- logN0_eig(
      z,
      U_curr,
      gam_curr,
      nugget = delta / 2,
      jitter = jitter
    )
    logZ_prop <- logN0_eig(
      z,
      U_prop,
      gam_prop,
      nugget = delta / 2,
      jitter = jitter
    )
    log_alpha <- (lf_prop - lf) +
      (gzy - gzx) +
      (logZ_prop - logZ_curr) +
      (lp_prop - lp_curr)
    acc <- log(runif(1)) < min(0, log_alpha)
    if (acc) {
      list(
        f = f_prop,
        lf = lf_prop,
        gf = gf_prop,
        theta = theta_prop,
        U = U_prop,
        gam = gam_prop,
        acc = TRUE
      )
    } else {
      list(
        f = f,
        lf = lf,
        gf = gf,
        theta = theta,
        U = U_curr,
        gam = gam_curr,
        acc = FALSE
      )
    }
  }
  
  tail_theta <- if (state_tail_length > 0L) {
    matrix(
      NA_real_,
      nrow = state_tail_length,
      ncol = 2L,
      dimnames = list(NULL, c("log_ell_cov", "log_sigma_cov"))
    )
  } else {
    NULL
  }
  tail_beta <- if (state_tail_length > 0L) {
    matrix(
      NA_real_,
      nrow = state_tail_length,
      ncol = p_coef,
      dimnames = list(NULL, colnames(Z_mean))
    )
  } else {
    NULL
  }
  tail_tau2 <- if (state_tail_length > 0L) {
    numeric(state_tail_length)
  } else {
    NULL
  }
  tail_aux_tau <- if (state_tail_length > 0L) {
    numeric(state_tail_length)
  } else {
    NULL
  }
  tail_eta <- if (state_tail_length > 0L) {
    numeric(state_tail_length)
  } else {
    NULL
  }
  tail_start_iteration <- if (state_tail_length > 0L) {
    n_burnin - state_tail_length + 1L
  } else {
    Inf
  }
  tail_index <- 0L
  
  store_tail_state <- function(global_warmup_iteration) {
    if (state_tail_length > 0L &&
        global_warmup_iteration >= tail_start_iteration) {
      tail_index <<- tail_index + 1L
      tail_theta[tail_index, ] <<- theta_cov
      tail_beta[tail_index, ] <<- coef_vec
      tail_tau2[tail_index] <<- tau2_cov
      tail_aux_tau[tail_index] <<- aux_tau_cov
      tail_eta[tail_index] <<- eta
    }
    invisible(NULL)
  }
  
  start_time <- Sys.time()
  if (verbose) {
    cat(sprintf(
      "\n[%s] %s exact sampler started.\n",
      timestamp_now(),
      phase_label
    ))
    cat("Observations          :", n, "\n")
    cat("Covariate columns     :", p_mean, "\n")
    cat("Blocked warm-up       :", blocked_iters, "\n")
    cat("Joint tuning          :", joint_tune_iters, "\n")
    cat("Posterior sampling    :", n_sample, "\n")
    cat("Thinning              :", thin, "\n")
    cat("Retained draws        :", n_saved_expected, "\n")
    cat("Saved pilot tail      :", state_tail_length, "\n")
  }
  
  ## Phase 1: blocked warm-up.
  acc_blocked_sum <- c(f = 0L, theta = 0L)
  eta_moved_blocked_sum <- 0L
  win_len <- 0L
  win_acc <- c(f = 0L, theta = 0L)
  win_id_blocked <- 0L
  theta_path_blocked <- matrix(
    NA_real_,
    nrow = max(blocked_iters, 1L),
    ncol = 2L,
    dimnames = list(NULL, c("log_ell_cov", "log_sigma_cov"))
  )
  cor_adapt_trace <- data.frame(
    iteration = integer(0),
    rho_cov = numeric(0)
  )
  
  if (blocked_iters > 0L) {
    for (i in seq_len(blocked_iters)) {
      out_f <- agrad_f_step(
        f,
        mu_vec,
        lf,
        gf,
        eta,
        U_sum,
        gam_sum,
        delta_state
      )
      f <- out_f$f
      lf <- out_f$lf
      gf <- out_f$gf
      acc_f <- as.integer(out_f$acc)
      
      out_theta <- theta_block_step(
        theta_cov,
        f,
        mu_vec,
        U_sum,
        gam_sum,
        kappa_cov,
        proposal_cor
      )
      theta_cov <- out_theta$theta
      U_sum <- out_theta$U
      gam_sum <- out_theta$gam
      acc_theta <- as.integer(out_theta$acc)
      
      upd <- update_eta_coef_tau(
        f = f,
        eta_current = eta,
        U = U_sum,
        gam = gam_sum,
        coef_vec = coef_vec,
        tau2_cov = tau2_cov,
        aux_tau_cov = aux_tau_cov
      )
      eta <- upd$eta
      lf <- upd$lf
      gf <- upd$gf
      coef_vec <- upd$coef_vec
      tau2_cov <- upd$tau2_cov
      aux_tau_cov <- upd$aux_tau_cov
      mu_vec <- upd$mu_vec
      
      acc_now <- c(f = acc_f, theta = acc_theta)
      acc_blocked_sum <- acc_blocked_sum + acc_now
      eta_moved_blocked_sum <-
        eta_moved_blocked_sum + upd$eta_moved
      theta_path_blocked[i, ] <- theta_cov
      
      if (adapt && adapt_block_cor && i >= cor_adapt_start &&
          i %% cor_adapt_interval == 0L) {
        idx_cor <- seq.int(max(1L, i - cor_adapt_window + 1L), i)
        proposal_cor <- estimate_cor_2x2(
          theta_path_blocked[idx_cor, , drop = FALSE],
          fallback = proposal_cor,
          shrinkage = cor_shrinkage,
          max_abs = cor_max_abs,
          min_n = min(200L, length(idx_cor))
        )
        cor_adapt_trace <- rbind(
          cor_adapt_trace,
          data.frame(iteration = i, rho_cov = proposal_cor[1L, 2L])
        )
      }
      
      if (adapt && i >= adapt_start) {
        win_len <- win_len + 1L
        win_acc <- win_acc + acc_now
        if (win_len == adapt_interval) {
          win_id_blocked <- win_id_blocked + 1L
          ar <- win_acc / win_len
          delta_state <- adapt_log_step(
            delta_state,
            ar["f"],
            target_f,
            win_id_blocked,
            adapt_rate_delta,
            t0_adapt,
            power_adapt,
            delta_min,
            delta_max
          )
          kappa_cov <- adapt_log_step(
            kappa_cov,
            ar["theta"],
            target_theta,
            win_id_blocked,
            adapt_rate_theta,
            t0_adapt,
            power_adapt,
            kappa_min,
            kappa_max
          )
          win_len <- 0L
          win_acc[] <- 0L
        }
      }
      
      store_tail_state(i)
      
      if (verbose && burnin_progress_every > 0L &&
          i %% burnin_progress_every == 0L) {
        progress_message(
          paste0(phase_label, " blocked warm-up"),
          i,
          blocked_iters
        )
      }
    }
  }
  
  acc_blocked <- acc_blocked_sum / max(1L, blocked_iters)
  eta_moved_blocked <-
    eta_moved_blocked_sum / max(1L, blocked_iters)
  proposal_cor_final <- sanitize_cor_2x2(
    proposal_cor,
    max_abs = cor_max_abs
  )
  
  theta_tail_blocked <- if (blocked_iters > 0L) {
    idx_start <- max(1L, blocked_iters - emp_cov_tail + 1L)
    theta_path_blocked[idx_start:blocked_iters, , drop = FALSE]
  } else {
    matrix(theta_cov, nrow = 1L)
  }
  
  if (!is.null(joint_cov_base_init)) {
    Sigma_theta_base <- sanitize_covariance_matrix(
      as.matrix(joint_cov_base_init) + cov_jitter * diag(2L)
    )
    if (!all(dim(Sigma_theta_base) == c(2L, 2L))) {
      stop("joint_cov_base_init must be 2 x 2.")
    }
    joint_cor_base <- sanitize_cor_2x2(
      safe_cor_from_cov(Sigma_theta_base),
      max_abs = cor_max_abs
    )
    joint_geometry_source <- "supplied_global_pilot_covariance"
  } else {
    joint_cor_base <- estimate_cor_2x2(
      theta_tail_blocked,
      fallback = proposal_cor_final,
      shrinkage = cor_shrinkage,
      max_abs = cor_max_abs,
      min_n = min(200L, nrow(theta_tail_blocked))
    )
    Sigma_theta_base <- proposal_cov_scaled_2x2(
      kappa_scalar = kappa_cov,
      prior_sd = prior_sd_theta,
      cor_mat = joint_cor_base,
      jitter = cov_jitter
    )
    Sigma_theta_base <- sanitize_covariance_matrix(Sigma_theta_base)
    joint_geometry_source <- "blocked_path_prior_scaled_correlation"
  }
  dimnames(joint_cor_base) <- list(
    names(prior_sd_theta),
    names(prior_sd_theta)
  )
  dimnames(Sigma_theta_base) <- list(
    names(prior_sd_theta),
    names(prior_sd_theta)
  )
  
  delta_joint <- delta_state
  
  if (verbose) {
    cat(sprintf(
      "[%s] %s blocked phase complete.\n",
      timestamp_now(),
      phase_label
    ))
    cat("  delta_joint          :", delta_joint, "\n")
    cat("  kappa_cov            :", kappa_cov, "\n")
    cat("  rho_cov              :", joint_cor_base[1L, 2L], "\n")
    cat("  joint geometry source:", joint_geometry_source, "\n")
  }
  
  ## Phase 2: joint tuning.
  acc_joint_tune_sum <- 0L
  acc_f_refresh_tune_sum <- 0L
  eta_moved_joint_tune_sum <- 0L
  joint_win_len <- 0L
  joint_win_acc <- 0L
  win_id_joint <- 0L
  
  if (joint_tune_iters > 0L) {
    for (i in seq_len(joint_tune_iters)) {
      Sigma_theta_eff <- (joint_theta_scale^2) * Sigma_theta_base
      out <- joint_step(
        f,
        mu_vec,
        lf,
        gf,
        theta_cov,
        eta,
        U_sum,
        gam_sum,
        delta_joint,
        Sigma_theta_eff
      )
      f <- out$f
      lf <- out$lf
      gf <- out$gf
      theta_cov <- out$theta
      U_sum <- out$U
      gam_sum <- out$gam
      acc_joint <- as.integer(out$acc)
      acc_joint_tune_sum <- acc_joint_tune_sum + acc_joint
      
      if (isTRUE(f_refresh_after_joint)) {
        out_f <- agrad_f_step(
          f,
          mu_vec,
          lf,
          gf,
          eta,
          U_sum,
          gam_sum,
          delta_joint
        )
        f <- out_f$f
        lf <- out_f$lf
        gf <- out_f$gf
        acc_f_refresh_tune_sum <- acc_f_refresh_tune_sum +
          as.integer(out_f$acc)
      }
      
      upd <- update_eta_coef_tau(
        f = f,
        eta_current = eta,
        U = U_sum,
        gam = gam_sum,
        coef_vec = coef_vec,
        tau2_cov = tau2_cov,
        aux_tau_cov = aux_tau_cov
      )
      eta <- upd$eta
      lf <- upd$lf
      gf <- upd$gf
      coef_vec <- upd$coef_vec
      tau2_cov <- upd$tau2_cov
      aux_tau_cov <- upd$aux_tau_cov
      mu_vec <- upd$mu_vec
      eta_moved_joint_tune_sum <-
        eta_moved_joint_tune_sum + upd$eta_moved
      
      if (adapt && i >= adapt_start) {
        joint_win_len <- joint_win_len + 1L
        joint_win_acc <- joint_win_acc + acc_joint
        if (joint_win_len == adapt_interval) {
          win_id_joint <- win_id_joint + 1L
          ar_joint <- joint_win_acc / joint_win_len
          joint_theta_scale <- adapt_log_step(
            joint_theta_scale,
            ar_joint,
            target_joint,
            win_id_joint,
            adapt_rate_joint_theta,
            t0_adapt,
            power_adapt,
            joint_theta_scale_min,
            joint_theta_scale_max
          )
          joint_win_len <- 0L
          joint_win_acc <- 0L
        }
      }
      
      store_tail_state(blocked_iters + i)
      
      if (verbose && burnin_progress_every > 0L &&
          i %% burnin_progress_every == 0L) {
        progress_message(
          paste0(phase_label, " joint tuning"),
          i,
          joint_tune_iters
        )
      }
    }
  }
  
  acc_joint_tune <- acc_joint_tune_sum / max(1L, joint_tune_iters)
  acc_f_refresh_tune <- acc_f_refresh_tune_sum / max(1L, joint_tune_iters)
  eta_moved_joint_tune <-
    eta_moved_joint_tune_sum / max(1L, joint_tune_iters)
  delta_final <- delta_joint
  joint_theta_scale_final <- joint_theta_scale
  Sigma_theta_final <- sanitize_covariance_matrix(
    (joint_theta_scale_final^2) * Sigma_theta_base
  )
  dimnames(Sigma_theta_final) <- list(
    names(prior_sd_theta),
    names(prior_sd_theta)
  )
  
  Sigma_block_final <- proposal_cov_scaled_2x2(
    kappa_scalar = kappa_cov,
    prior_sd = prior_sd_theta,
    cor_mat = proposal_cor_final,
    jitter = cov_jitter
  )
  Sigma_block_final <- sanitize_covariance_matrix(Sigma_block_final)
  dimnames(Sigma_block_final) <- list(
    names(prior_sd_theta),
    names(prior_sd_theta)
  )
  
  if (state_tail_length > 0L && tail_index != state_tail_length) {
    stop(
      "The pilot tail was not filled exactly: got ",
      tail_index,
      ", expected ",
      state_tail_length,
      "."
    )
  }
  
  ## Phase 3: posterior sampling with frozen proposal geometry.
  score_sum <- numeric(n)
  n_saved <- 0L
  
  predictive_T_draws <- if (n_saved_expected > 0L) {
    matrix(
      NA_integer_,
      nrow = n_saved_expected,
      ncol = n
    )
  } else {
    matrix(integer(0), nrow = 0L, ncol = n)
  }
  
  if (n_saved_expected > 0L) {
    if (
      is.null(predictive_seed) ||
      length(predictive_seed) != 1L ||
      !is.finite(predictive_seed)
    ) {
      stop(
        "predictive_seed is required for GP-Cov + PU predictive labels."
      )
    }
    
    largest_predictive_seed <- as.numeric(predictive_seed) +
      n_saved_expected
    
    if (
      predictive_seed < 1 ||
      largest_predictive_seed > .Machine$integer.max
    ) {
      stop("The predictive-label seed stream exceeds the R integer range.")
    }
  }
  
  eta_draws <- numeric(n_saved_expected)
  max_probability_draws <- numeric(n_saved_expected)
  q99_probability_draws <- numeric(n_saved_expected)
  top5_mean_probability_draws <- numeric(n_saved_expected)
  beta0_draws <- numeric(n_saved_expected)
  tau_cov_draws <- numeric(n_saved_expected)
  theta_draws <- if (n_saved_expected > 0L) {
    matrix(
      NA_real_,
      nrow = n_saved_expected,
      ncol = 2L,
      dimnames = list(NULL, c("ell_cov", "sigma_cov"))
    )
  } else {
    matrix(numeric(0), nrow = 0L, ncol = 2L)
  }
  
  acc_joint_sampling_sum <- 0L
  acc_f_refresh_sampling_sum <- 0L
  eta_moved_sampling_sum <- 0L
  
  if (n_sample > 0L) {
    for (t in seq_len(n_sample)) {
      out <- joint_step(
        f,
        mu_vec,
        lf,
        gf,
        theta_cov,
        eta,
        U_sum,
        gam_sum,
        delta_final,
        Sigma_theta_final
      )
      f <- out$f
      lf <- out$lf
      gf <- out$gf
      theta_cov <- out$theta
      U_sum <- out$U
      gam_sum <- out$gam
      acc_joint_sampling_sum <- acc_joint_sampling_sum +
        as.integer(out$acc)
      
      if (isTRUE(f_refresh_after_joint)) {
        out_f <- agrad_f_step(
          f,
          mu_vec,
          lf,
          gf,
          eta,
          U_sum,
          gam_sum,
          delta_final
        )
        f <- out_f$f
        lf <- out_f$lf
        gf <- out_f$gf
        acc_f_refresh_sampling_sum <- acc_f_refresh_sampling_sum +
          as.integer(out_f$acc)
      }
      
      upd <- update_eta_coef_tau(
        f = f,
        eta_current = eta,
        U = U_sum,
        gam = gam_sum,
        coef_vec = coef_vec,
        tau2_cov = tau2_cov,
        aux_tau_cov = aux_tau_cov
      )
      eta <- upd$eta
      lf <- upd$lf
      gf <- upd$gf
      coef_vec <- upd$coef_vec
      tau2_cov <- upd$tau2_cov
      aux_tau_cov <- upd$aux_tau_cov
      mu_vec <- upd$mu_vec
      eta_moved_sampling_sum <-
        eta_moved_sampling_sum + upd$eta_moved
      
      if (t %% thin == 0L) {
        n_saved <- n_saved + 1L
        p_current <- plogis(f)
        score_sum <- score_sum + p_current
        
        ## For every unit and retained posterior draw, construct the
        ## predictive latent label from T_i ~ Bernoulli(plogis(f_i)).
        predictive_draw_seed <- as.integer(
          as.numeric(predictive_seed) + n_saved
        )
        
        predictive_T_draws[n_saved, ] <- with_local_seed(
          predictive_draw_seed,
          {
            stats::rbinom(
              n,
              size = 1L,
              prob = p_current
            )
          }
        )
        
        eta_draws[n_saved] <- eta
        max_probability_draws[n_saved] <- max(p_current)
        q99_probability_draws[n_saved] <- unname(
          quantile(p_current, 0.99, names = FALSE)
        )
        top5_mean_probability_draws[n_saved] <- mean(
          sort(p_current, decreasing = TRUE)[1:5]
        )
        beta0_draws[n_saved] <- coef_vec[1L]
        tau_cov_draws[n_saved] <- sqrt(tau2_cov)
        theta_draws[n_saved, ] <- exp(theta_cov)
      }
      
      if (verbose && sample_progress_every > 0L &&
          t %% sample_progress_every == 0L) {
        progress_message(
          paste0(phase_label, " posterior sampling"),
          t,
          n_sample
        )
      }
    }
  }
  
  if (n_saved != n_saved_expected) {
    stop(
      "Retained-draw count mismatch: got ",
      n_saved,
      ", expected ",
      n_saved_expected,
      "."
    )
  }
  
  end_time <- Sys.time()
  runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
  acc_sampling <- acc_joint_sampling_sum / max(1L, n_sample)
  acc_f_refresh_sampling <- acc_f_refresh_sampling_sum / max(1L, n_sample)
  eta_moved_sampling <-
    eta_moved_sampling_sum / max(1L, n_sample)
  
  if (verbose) {
    cat(sprintf(
      "[%s] %s exact sampler finished.\n",
      timestamp_now(),
      phase_label
    ))
    cat("Runtime                  :", format_duration(runtime_sec), "\n")
    cat("Blocked acceptance f     :", acc_blocked["f"], "\n")
    cat("Blocked acceptance theta :", acc_blocked["theta"], "\n")
    cat("Joint tune acceptance    :", acc_joint_tune, "\n")
    cat("Joint sample acceptance  :", acc_sampling, "\n")
    cat("f-refresh sample accept. :", acc_f_refresh_sampling, "\n")
    cat("eta move rate            :", eta_moved_sampling, "\n")
    cat("Final ell / sigma        :", exp(theta_cov[1L]), "/",
        exp(theta_cov[2L]), "\n")
    cat("Final eta                :", eta, "\n")
    cat("Final tau_cov            :", sqrt(tau2_cov), "\n")
    cat("Final delta              :", delta_final, "\n")
    cat("Final joint theta scale  :", joint_theta_scale_final, "\n")
  }
  
  list(
    state = list(
      theta_log = setNames(
        as.numeric(theta_cov),
        c("log_ell_cov", "log_sigma_cov")
      ),
      beta = setNames(as.numeric(coef_vec), colnames(Z_mean)),
      tau2_cov = as.numeric(tau2_cov),
      aux_tau_cov = as.numeric(aux_tau_cov),
      eta = as.numeric(eta),
      f = as.numeric(f)
    ),
    tail = if (state_tail_length > 0L) {
      list(
        theta_log = tail_theta,
        beta = tail_beta,
        tau2_cov = tail_tau2,
        aux_tau_cov = tail_aux_tau,
        eta = tail_eta
      )
    } else {
      NULL
    },
    posterior_mean_scores = if (n_saved > 0L) score_sum / n_saved else NULL,
    predictive_T_draws = predictive_T_draws,
    eta_draws = eta_draws,
    max_probability_draws = max_probability_draws,
    q99_probability_draws = q99_probability_draws,
    top5_mean_probability_draws = top5_mean_probability_draws,
    beta0_draws = beta0_draws,
    tau_cov_draws = tau_cov_draws,
    theta_draws = theta_draws,
    n_saved = n_saved,
    proposal = list(
      delta_final = as.numeric(delta_final),
      kappa_cov_final = as.numeric(kappa_cov),
      proposal_cor_final = proposal_cor_final,
      block_cov_final = Sigma_block_final,
      joint_cor_base = joint_cor_base,
      Sigma_theta_base = Sigma_theta_base,
      joint_theta_scale_final = as.numeric(joint_theta_scale_final),
      Sigma_theta_final = Sigma_theta_final,
      joint_geometry_source = joint_geometry_source,
      cor_adapt_trace = cor_adapt_trace
    ),
    acceptance = list(
      blocked_f = unname(acc_blocked["f"]),
      blocked_theta = unname(acc_blocked["theta"]),
      joint_tune = acc_joint_tune,
      f_refresh_tune = acc_f_refresh_tune,
      joint_sampling = acc_sampling,
      f_refresh_sampling = acc_f_refresh_sampling,
      eta_moved_blocked = eta_moved_blocked,
      eta_moved_joint_tune = eta_moved_joint_tune,
      eta_moved_sampling = eta_moved_sampling
    ),
    phases = c(
      blocked_warmup = blocked_iters,
      joint_tune = joint_tune_iters,
      posterior_sampling = n_sample,
      retained = n_saved
    ),
    time_start = start_time,
    time_end = end_time,
    time_fit_sec = runtime_sec
  )
}

## ============================================================================

## ============================================================================
## 9. COVARIATE, SAR-LABEL AND EXACT-PRIOR BUNDLES
## ============================================================================

validate_sar_bundle_labels <- function(
    observed_y,
    hidden_positive_idx,
    observed_positive_idx,
    signal_type,
    context_label
) {
  observed_y <- as.integer(observed_y)
  hidden_positive_idx <- sort(unique(as.integer(hidden_positive_idx)))
  observed_positive_idx <- sort(unique(as.integer(observed_positive_idx)))

  if (
    length(observed_y) != N_NODES ||
    !all(observed_y %in% c(0L, 1L)) ||
    sum(observed_y) != N_OBSERVED_POS
  ) {
    stop(context_label, " has an invalid observed-label vector.")
  }

  if (
    length(hidden_positive_idx) != N_HIDDEN ||
    length(observed_positive_idx) != N_OBSERVED_POS ||
    any(!hidden_positive_idx %in% TRUE_POS_IDX) ||
    any(!observed_positive_idx %in% TRUE_POS_IDX) ||
    !identical(
      sort(c(hidden_positive_idx, observed_positive_idx)),
      TRUE_POS_IDX
    ) ||
    any(observed_y[hidden_positive_idx] != 0L) ||
    any(observed_y[observed_positive_idx] != 1L)
  ) {
    stop(context_label, " has an inconsistent SAR positive-label partition.")
  }

  if (
    length(signal_type) != N_NODES ||
    any(is.na(signal_type))
  ) {
    stop(context_label, " has an invalid signal_type vector.")
  }

  invisible(TRUE)
}

finalize_covariate_bundle <- function(
    source_label,
    X_cov,
    X_raw,
    seeds,
    observed_y,
    hidden_positive_idx,
    observed_positive_idx,
    signal_type,
    sar_selection_summary,
    sar_selection_diagnostics,
    oracle_covariate_score = NULL,
    oracle_covariate_block_scores = NULL,
    standardized_positive_oracle_score = NULL,
    positive_hidden_weights = NULL,
    covariate_profile = NULL,
    X_center = NULL,
    X_scale = NULL,
    calibration_label
) {
  X_cov <- as.matrix(X_cov)
  X_raw <- as.matrix(X_raw)

  if (
    nrow(X_cov) != N_NODES ||
    ncol(X_cov) != COVARIATE_MEAN_DIMENSION ||
    any(!is.finite(X_cov))
  ) {
    stop("A covariate bundle has incompatible standardized covariates.")
  }

  if (
    nrow(X_raw) != N_NODES ||
    ncol(X_raw) != COVARIATE_MEAN_DIMENSION ||
    any(!is.finite(X_raw))
  ) {
    stop("A covariate bundle has incompatible raw covariates.")
  }

  validate_sar_bundle_labels(
    observed_y = observed_y,
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = observed_positive_idx,
    signal_type = signal_type,
    context_label = source_label
  )

  if (is.null(colnames(X_cov))) {
    colnames(X_cov) <- paste0(
      "X",
      seq_len(ncol(X_cov))
    )
  }

  colnames(X_cov) <- make.unique(
    colnames(X_cov)
  )

  prior_calibration <- calibrate_covariate_gp_prior(
    X_cov = X_cov,
    calibration_label = calibration_label
  )

  list(
    source = source_label,
    X_cov = X_cov,
    X_raw = X_raw,
    seeds = seeds,
    observed_y = as.integer(observed_y),
    hidden_positive_idx = sort(as.integer(hidden_positive_idx)),
    observed_positive_idx = sort(as.integer(observed_positive_idx)),
    signal_type = as.character(signal_type),
    sar_selection_summary = sar_selection_summary,
    sar_selection_diagnostics = sar_selection_diagnostics,
    oracle_covariate_score = oracle_covariate_score,
    oracle_covariate_block_scores = oracle_covariate_block_scores,
    standardized_positive_oracle_score =
      standardized_positive_oracle_score,
    positive_hidden_weights = positive_hidden_weights,
    covariate_profile = covariate_profile,
    X_center = X_center,
    X_scale = X_scale,
    prior_calibration = prior_calibration
  )
}

make_generated_covariate_bundle <- function(
    seed_row,
    source_label,
    calibration_label
) {
  required_seed_columns <- c(
    "covariate_seed",
    "sar_seed"
  )

  if (!all(required_seed_columns %in% names(seed_row))) {
    stop(
      "The generated SAR covariate bundle requires seeds: ",
      paste(required_seed_columns, collapse = ", "),
      "."
    )
  }

  covariate_set <- make_covariate_set(
    seed = seed_row$covariate_seed,
    signal_idx = COVARIATE_SIGNAL_IDX
  )

  sar_design <- build_sar_label_design(
    X_raw = covariate_set$X_raw,
    sar_seed = seed_row$sar_seed
  )

  finalize_covariate_bundle(
    source_label = source_label,
    X_cov = covariate_set$X_std,
    X_raw = covariate_set$X_raw,
    seeds = seed_row,
    observed_y = sar_design$observed_y,
    hidden_positive_idx = sar_design$hidden_positive_idx,
    observed_positive_idx = sar_design$observed_positive_idx,
    signal_type = sar_design$signal_type,
    sar_selection_summary = sar_design$selection_summary,
    sar_selection_diagnostics = sar_design$selection_diagnostics,
    oracle_covariate_score = sar_design$oracle_score,
    oracle_covariate_block_scores = sar_design$oracle_block_scores,
    standardized_positive_oracle_score =
      sar_design$standardized_positive_score,
    positive_hidden_weights = sar_design$hidden_weights,
    covariate_profile = covariate_set$covariate_profile,
    X_center = covariate_set$center,
    X_scale = covariate_set$scale,
    calibration_label = calibration_label
  )
}

load_or_generate_production_bundle <- function(
    replicate_number
) {
  seed_row <- production_seed_table[
    production_seed_table$replicate ==
      replicate_number,
    ,
    drop = FALSE
  ]

  if (nrow(seed_row) != 1L) {
    stop(
      "Production seed lookup failed for replicate ",
      replicate_number,
      "."
    )
  }

  calibration_label <- sprintf(
    "production_%03d",
    replicate_number
  )

  geometry_file <- file.path(
    COMPETING_GEOMETRY_DIR,
    sprintf(
      "replicate_%03d_geometry.rds",
      replicate_number
    )
  )

  if (
    isTRUE(USE_SAVED_COMPETING_GEOMETRIES) &&
    file.exists(geometry_file)
  ) {
    geometry <- readRDS(geometry_file)

    required_fields <- c(
      "T",
      "Y",
      "X_raw",
      "X_std",
      "hidden_positive_idx",
      "observed_positive_idx",
      "signal_type",
      "sar_selection_summary",
      "sar_selection_diagnostics"
    )

    missing_fields <- required_fields[
      !vapply(
        required_fields,
        function(name_now) {
          !is.null(geometry[[name_now]])
        },
        logical(1)
      )
    ]

    if (length(missing_fields) > 0L) {
      stop(
        "Saved competing SAR geometry is missing: ",
        paste(missing_fields, collapse = ", ")
      )
    }

    if (!identical(
      as.integer(geometry$T),
      TRUE_T
    )) {
      stop(
        "Saved competing SAR geometry has incompatible latent truth T."
      )
    }

    if (!is.null(geometry$seeds)) {
      for (seed_name in c("covariate_seed", "sar_seed")) {
        if (
          seed_name %in% names(geometry$seeds) &&
          !identical(
            as.integer(geometry$seeds[[seed_name]][1L]),
            as.integer(seed_row[[seed_name]])
          )
        ) {
          stop(
            "Saved competing SAR geometry has an incompatible ",
            seed_name,
            "."
          )
        }
      }
    }

    return(finalize_covariate_bundle(
      source_label =
        "saved_competing_SAR_geometry",
      X_cov = geometry$X_std,
      X_raw = geometry$X_raw,
      seeds = seed_row,
      observed_y = geometry$Y,
      hidden_positive_idx =
        geometry$hidden_positive_idx,
      observed_positive_idx =
        geometry$observed_positive_idx,
      signal_type = geometry$signal_type,
      sar_selection_summary =
        geometry$sar_selection_summary,
      sar_selection_diagnostics =
        geometry$sar_selection_diagnostics,
      oracle_covariate_score =
        geometry$oracle_covariate_score,
      oracle_covariate_block_scores =
        geometry$oracle_covariate_block_scores,
      standardized_positive_oracle_score =
        geometry$standardized_positive_oracle_score,
      positive_hidden_weights =
        geometry$positive_hidden_weights,
      covariate_profile =
        geometry$covariate_profile,
      X_center = geometry$X_center,
      X_scale = geometry$X_scale,
      calibration_label =
        calibration_label
    ))
  }

  make_generated_covariate_bundle(
    seed_row = seed_row,
    source_label =
      "regenerated_from_competing_SAR_seeds",
    calibration_label =
      calibration_label
  )
}

## ============================================================================
## 10. FIVE INDEPENDENT GLOBAL PILOTS
## ============================================================================

PRIOR_CALIBRATION_RULES <- list(
  correlation_mode =
    "relative_to_current_attainable_range",
  corr_pct_range_all =
    CORR_PCT_RANGE_ALL,
  ell_interval_probs =
    ELL_INTERVAL_PROBS,
  sigma_interval_probs =
    SIGMA_INTERVAL_PROBS,
  effective_sd_center_cov =
    EFFECTIVE_SD_CENTER_COV,
  sigma_range_factor =
    SIGMA_RANGE_FACTOR,
  ell_bounds_cov =
    ELL_BOUNDS_COV,
  n_ell_grid_coarse =
    N_ELL_GRID_COARSE,
  root_tol =
    ROOT_TOL,
  half_t_df =
    HALFT_DF,
  mean_rms_95_mult =
    MEAN_RMS_95_MULT,
  covariate_mean_dimension =
    COVARIATE_MEAN_DIMENSION,
  b_cov_fixed =
    B_COV_FIXED,
  tau_scale_cov_unrounded =
    TAU_SCALE_COV_UNROUNDED,
  tau_scale_cov_active =
    TAU_SCALE_COV_ACTIVE,
  round_active_priors =
    ROUND_ACTIVE_PRIORS,
  prior_round_digits =
    PRIOR_ROUND_DIGITS,
  beta0_prior =
    beta0_prior_for_sampler,
  eta_prior =
    eta_prior_for_sampler
)

PILOT_SETTINGS_SIGNATURE <- list(
  script_version = SCRIPT_VERSION,
  master_seed = MASTER_SEED,
  n_pilot = N_PILOT,
  pilot_burnin = PILOT_BURNIN,
  pilot_tail_length = PILOT_TAIL_LENGTH,
  pilot_blocked_frac = PILOT_BLOCKED_FRAC,

  data_design = list(
    true_t = TRUE_T,
    true_zero_idx = TRUE_ZERO_IDX,
    true_positive_idx = TRUE_POS_IDX,
    covariates_only_idx = COV_ONLY_IDX,
    network_only_idx = NET_ONLY_IDX,
    both_idx = BOTH_IDX,
    covariate_signal_idx =
      COVARIATE_SIGNAL_IDX,
    network_signal_idx =
      NETWORK_SIGNAL_IDX,
    n_hidden_per_pilot =
      N_HIDDEN,
    n_observed_positive_per_pilot =
      N_OBSERVED_POS
  ),

  label_mechanism = list(
    type = "covariate-driven SAR",
    exactly_n_hidden_without_replacement =
      N_HIDDEN,
    lambda = LAMBDA_SAR,
    positive_score_standardization =
      "within the 45 true positives",
    z_clip = SAR_Z_CLIP,
    log_hiding_weight =
      "-lambda * clipped standardized oracle covariate score",
    oracle_score_available_to_fitted_model =
      FALSE
  ),

  covariate_dgp = list(
    mu_linear = MU_LINEAR,
    sd_linear = SD_LINEAR,
    block2_x3_offset =
      BLOCK2_X3_OFFSET,
    block2_x4_center =
      BLOCK2_X4_CENTER,
    block2_x4_offset =
      BLOCK2_X4_OFFSET,
    block2_zero_sd_x3 =
      BLOCK2_ZERO_SD_X3,
    block2_zero_sd_x4 =
      BLOCK2_ZERO_SD_X4,
    block2_pos_sd_x3 =
      BLOCK2_POS_SD_X3,
    block2_pos_sd_x4 =
      BLOCK2_POS_SD_X4,
    disc_radius = DISC_RADIUS,
    ring_radius = RING_RADIUS,
    ring_sd = RING_SD,
    parabola_coef = PARABOLA_COEF,
    delta_parabola = DELTA_PARABOLA,
    sd_parabola_zero =
      SD_PARABOLA_ZERO,
    sd_parabola_pos =
      SD_PARABOLA_POS,
    rho_zero = RHO_ZERO,
    rho_pos = RHO_POS
  ),

  prior_calibration_rules =
    PRIOR_CALIBRATION_RULES,

  sampler_controls = list(
    min_blocked_warmup =
      MIN_BLOCKED_WARMUP,
    max_log_sigma_cov =
      MAX_LOG_SIGMA_COV,
    min_log_sigma =
      MIN_LOG_SIGMA,
    logl_cov_bounds =
      LOGL_COV_BOUNDS,
    eta_slice_width =
      ETA_SLICE_WIDTH,
    targets = c(
      f = TARGET_F,
      theta = TARGET_THETA,
      joint = TARGET_JOINT
    )
  )
)

run_or_load_one_pilot <- function(
    pilot_number
) {
  checkpoint_file <- file.path(
    PILOT_CHECKPOINT_DIR,
    sprintf(
      "GP_Cov_PU_pilot_%02d.rds",
      pilot_number
    )
  )
  
  if (isTRUE(RESUME_PILOTS) &&
      !isTRUE(FORCE_FRESH_PILOTS) &&
      file.exists(checkpoint_file)) {
    candidate <- tryCatch(
      readRDS(checkpoint_file),
      error = function(e) NULL
    )
    
    if (is.list(candidate) &&
        identical(
          candidate$settings_signature,
          PILOT_SETTINGS_SIGNATURE
        ) &&
        identical(
          as.integer(candidate$pilot),
          as.integer(pilot_number)
        ) &&
        isTRUE(candidate$complete)) {
      cat(sprintf(
        paste0(
          "[%s] Reusing self-calibrated ",
          "GP-Cov + PU SAR pilot %02d / %02d\n"
        ),
        timestamp_now(),
        pilot_number,
        N_PILOT
      ))
      
      return(candidate)
    }
  }
  
  seed_row <- pilot_seed_table[
    pilot_seed_table$pilot ==
      pilot_number,
    ,
    drop = FALSE
  ]
  
  if (nrow(seed_row) != 1L) {
    stop("Pilot seed lookup failed.")
  }
  
  bundle <- make_generated_covariate_bundle(
    seed_row = seed_row,
    source_label =
      "independent_global_SAR_pilot",
    calibration_label = sprintf(
      "pilot_%02d",
      pilot_number
    )
  )

  validate_sar_bundle_labels(
    observed_y = bundle$observed_y,
    hidden_positive_idx =
      bundle$hidden_positive_idx,
    observed_positive_idx =
      bundle$observed_positive_idx,
    signal_type = bundle$signal_type,
    context_label = paste0(
      "GP-Cov + PU SAR pilot ",
      pilot_number
    )
  )

  pilot_sar_summary <-
    bundle$sar_selection_summary

  pilot_sar_summary$pilot <-
    pilot_number

  pilot_sar_summary <-
    pilot_sar_summary[
      ,
      c(
        "scenario",
        "pilot",
        setdiff(
          names(pilot_sar_summary),
          c("scenario", "pilot")
        )
      ),
      drop = FALSE
    ]
  
  pilot_prior_summary <-
    bundle$prior_calibration$summary
  
  pilot_prior_summary$stage <- "pilot"
  pilot_prior_summary$index <- pilot_number
  
  pilot_prior_summary <- pilot_prior_summary[
    ,
    c(
      "stage",
      "index",
      setdiff(
        names(pilot_prior_summary),
        c("stage", "index")
      )
    ),
    drop = FALSE
  ]
  
  write.csv(
    pilot_prior_summary,
    file.path(
      PRIOR_CALIBRATION_DIR,
      sprintf(
        "GP_Cov_PU_pilot_%02d_prior_calibration.csv",
        pilot_number
      )
    ),
    row.names = FALSE
  )
  
  if (isTRUE(SAVE_PILOT_GEOMETRIES)) {
    saveRDS(
      list(
        settings_signature =
          PILOT_SETTINGS_SIGNATURE,
        pilot = pilot_number,
        seeds = seed_row,
        T = TRUE_T,
        Y = bundle$observed_y,
        hidden_positive_idx =
          bundle$hidden_positive_idx,
        observed_positive_idx =
          bundle$observed_positive_idx,
        positive_signal_group =
          POSITIVE_SIGNAL_GROUP,
        signal_type =
          bundle$signal_type,
        X_raw = bundle$X_raw,
        X_cov = bundle$X_cov,
        covariate_profile =
          bundle$covariate_profile,
        X_center = bundle$X_center,
        X_scale = bundle$X_scale,
        oracle_covariate_score =
          bundle$oracle_covariate_score,
        oracle_covariate_block_scores =
          bundle$oracle_covariate_block_scores,
        standardized_positive_oracle_score =
          bundle$standardized_positive_oracle_score,
        positive_hidden_weights =
          bundle$positive_hidden_weights,
        sar_selection_summary =
          pilot_sar_summary,
        sar_selection_diagnostics =
          bundle$sar_selection_diagnostics,
        prior_calibration =
          bundle$prior_calibration
      ),
      file.path(
        PILOT_GEOMETRY_DIR,
        sprintf(
          "GP_Cov_PU_pilot_%02d_geometry.rds",
          pilot_number
        )
      )
    )
  }
  
  mcmc_seed <- as.integer(
    PILOT_MCMC_SEED_BASE +
      10000L * pilot_number
  )
  
  cat("\n\n####################################################################\n")
  cat(sprintf(
    "# SELF-CALIBRATED GP-COV + PU SAR GLOBAL PILOT %d / %d\n",
    pilot_number,
    N_PILOT
  ))
  cat("####################################################################\n")
  cat("Covariate seed      :",
      seed_row$covariate_seed, "\n")
  cat("SAR seed            :",
      seed_row$sar_seed, "\n")
  cat(
    "SAR hidden composition: cov=",
    pilot_sar_summary$hidden_covariates_only,
    ", net=",
    pilot_sar_summary$hidden_network_only,
    ", both=",
    pilot_sar_summary$hidden_both,
    "\n",
    sep = ""
  )
  cat("MCMC seed           :",
      mcmc_seed, "\n")
  cat("Active kernel prior :\n")
  print(
    bundle$prior_calibration$active_hyperprior
  )
  cat("Active tau scale    :",
      bundle$prior_calibration$tau_scale,
      "\n")
  
  set.seed(mcmc_seed)
  
  fit <- run_gp_cov_pu_exact_chain(
    Y = bundle$observed_y,
    X_cov = bundle$X_cov,
    hyperprior_cov =
      bundle$prior_calibration$active_hyperprior,
    m_beta =
      beta0_prior_for_sampler$mean,
    s_beta =
      beta0_prior_for_sampler$sd,
    tau2_default =
      bundle$prior_calibration$tau2_default,
    scale_tau =
      bundle$prior_calibration$tau_scale,
    a_eta =
      eta_prior_for_sampler$a_eta,
    b_eta =
      eta_prior_for_sampler$b_eta,
    
    n_burnin = PILOT_BURNIN,
    n_sample = 0L,
    thin = 1L,
    blocked_frac = PILOT_BLOCKED_FRAC,
    min_blocked_warmup =
      MIN_BLOCKED_WARMUP,
    
    init_state = NULL,
    delta_init = 3.5,
    kappa_cov_init = KAPPA_COV_INIT,
    proposal_cor_init =
      PROPOSAL_COR_INIT,
    
    adapt_block_cor =
      ADAPT_BLOCK_CORRELATION_PILOTS,
    cor_adapt_start =
      COR_ADAPT_START,
    cor_adapt_interval =
      COR_ADAPT_INTERVAL,
    cor_adapt_window =
      COR_ADAPT_WINDOW,
    cor_shrinkage =
      COR_SHRINKAGE,
    cor_max_abs =
      COR_MAX_ABS,
    
    joint_cov_base_init = NULL,
    joint_theta_scale_init =
      JOINT_THETA_SCALE_INIT_PILOT,
    f_refresh_after_joint =
      F_REFRESH_AFTER_JOINT,
    eta_slice_width =
      ETA_SLICE_WIDTH,
    initial_eta =
      INITIAL_ETA,
    
    nu_tau = NU_TAU,
    sample_tau = SAMPLE_TAU,
    
    adapt = ADAPT,
    adapt_start = ADAPT_START,
    adapt_interval = ADAPT_INTERVAL,
    target_f = TARGET_F,
    target_theta = TARGET_THETA,
    target_joint = TARGET_JOINT,
    adapt_rate_delta =
      ADAPT_RATE_DELTA,
    adapt_rate_theta =
      ADAPT_RATE_THETA,
    adapt_rate_joint_theta =
      ADAPT_RATE_JOINT_THETA,
    t0_adapt = T0_ADAPT,
    power_adapt = POWER_ADAPT,
    
    delta_min = DELTA_MIN,
    delta_max = DELTA_MAX,
    kappa_min = KAPPA_MIN,
    kappa_max = KAPPA_MAX,
    joint_theta_scale_min =
      JOINT_THETA_SCALE_MIN,
    joint_theta_scale_max =
      JOINT_THETA_SCALE_MAX,
    
    emp_cov_tail = EMP_COV_TAIL,
    cov_jitter = COV_JITTER,
    max_log_sigma_cov =
      MAX_LOG_SIGMA_COV,
    min_log_sigma =
      MIN_LOG_SIGMA,
    logl_cov_bounds =
      LOGL_COV_BOUNDS,
    jitter = JITTER,
    
    state_tail_length =
      PILOT_TAIL_LENGTH,
    predictive_seed = NULL,
    
    phase_label = paste0(
      "GP-Cov + PU SAR pilot ",
      sprintf("%02d", pilot_number)
    ),
    burnin_progress_every =
      PILOT_PROGRESS_EVERY,
    sample_progress_every = 0L,
    verbose = TRUE
  )
  
  result <- list(
    settings_signature =
      PILOT_SETTINGS_SIGNATURE,
    pilot = pilot_number,
    complete = TRUE,
    geometry_source = bundle$source,
    covariate_seed =
      seed_row$covariate_seed,
    sar_seed =
      seed_row$sar_seed,
    observed_y =
      bundle$observed_y,
    hidden_positive_idx =
      bundle$hidden_positive_idx,
    observed_positive_idx =
      bundle$observed_positive_idx,
    sar_selection_summary =
      pilot_sar_summary,
    mcmc_seed = mcmc_seed,
    
    active_priors = list(
      hyperprior =
        bundle$prior_calibration$active_hyperprior,
      hyperprior_unrounded =
        bundle$prior_calibration$hyperprior_unrounded,
      tau_scale =
        bundle$prior_calibration$tau_scale,
      tau_scale_unrounded =
        bundle$prior_calibration$tau_scale_unrounded,
      tau2_default =
        bundle$prior_calibration$tau2_default,
      theta_prior_mean =
        bundle$prior_calibration$theta_prior_mean,
      theta_prior_sd =
        bundle$prior_calibration$theta_prior_sd
    ),
    
    prior_calibration_summary =
      pilot_prior_summary,
    tail = fit$tail,
    final_state = fit$state,
    final_proposal = fit$proposal,
    acceptance = fit$acceptance,
    runtime_sec = fit$time_fit_sec
  )
  
  saveRDS(
    result,
    checkpoint_file
  )
  
  rm(bundle, fit)
  invisible(gc())
  
  result
}

run_all_pilots <- function() {
  pilot_results <- vector(
    "list",
    N_PILOT
  )
  
  for (pilot_number in seq_len(N_PILOT)) {
    pilot_results[[pilot_number]] <-
      run_or_load_one_pilot(
        pilot_number
      )
  }
  
  pilot_results
}

build_global_calibration <- function(
    pilot_results
) {
  expected_rows <-
    N_PILOT * PILOT_TAIL_LENGTH
  
  theta_names <- c(
    "log_ell_cov",
    "log_sigma_cov"
  )
  
  pooled_theta_raw <- do.call(
    rbind,
    lapply(
      pilot_results,
      function(x) {
        x$tail$theta_log
      }
    )
  )
  
  pooled_theta_standardized <- do.call(
    rbind,
    lapply(
      pilot_results,
      function(x) {
        centered <- sweep(
          x$tail$theta_log,
          2L,
          x$active_priors$theta_prior_mean[
            theta_names
          ],
          "-"
        )
        
        sweep(
          centered,
          2L,
          x$active_priors$theta_prior_sd[
            theta_names
          ],
          "/"
        )
      }
    )
  )
  
  colnames(pooled_theta_standardized) <-
    theta_names
  
  pooled_beta <- do.call(
    rbind,
    lapply(
      pilot_results,
      function(x) {
        x$tail$beta
      }
    )
  )
  
  pooled_tau2_ratio <- unlist(
    lapply(
      pilot_results,
      function(x) {
        x$tail$tau2_cov /
          x$active_priors$tau_scale^2
      }
    ),
    use.names = FALSE
  )
  
  pooled_aux_ratio <- unlist(
    lapply(
      pilot_results,
      function(x) {
        x$tail$aux_tau_cov *
          x$active_priors$tau_scale^2
      }
    ),
    use.names = FALSE
  )
  
  pooled_eta <- unlist(
    lapply(
      pilot_results,
      function(x) {
        x$tail$eta
      }
    ),
    use.names = FALSE
  )
  
  if (nrow(pooled_theta_raw) != expected_rows ||
      nrow(pooled_theta_standardized) !=
      expected_rows ||
      nrow(pooled_beta) != expected_rows ||
      length(pooled_tau2_ratio) !=
      expected_rows ||
      length(pooled_aux_ratio) !=
      expected_rows ||
      length(pooled_eta) !=
      expected_rows) {
    stop(
      "The pooled pilot tails have an unexpected number of rows."
    )
  }
  
  pooled_theta_correlation <-
    estimate_cor_2x2(
      pooled_theta_standardized,
      fallback = diag(2),
      shrinkage =
        GLOBAL_COR_SHRINKAGE,
      max_abs = COR_MAX_ABS,
      min_n = 200L
    )
  
  standardized_final_joint_covariances <-
    lapply(
      pilot_results,
      function(x) {
        prior_sd <-
          x$active_priors$theta_prior_sd[
            theta_names
          ]
        
        D_inverse <- diag(
          1 / prior_sd,
          2L
        )
        
        sanitize_covariance_matrix(
          D_inverse %*%
            x$final_proposal$Sigma_theta_final %*%
            D_inverse
        )
      }
    )
  
  global_joint_covariance_standardized <-
    median_covariance_2x2(
      standardized_final_joint_covariances
    )
  
  dimnames(
    global_joint_covariance_standardized
  ) <- list(
    theta_names,
    theta_names
  )
  
  global_kappa <- median(vapply(
    pilot_results,
    function(x) {
      x$final_proposal$kappa_cov_final
    },
    numeric(1)
  ))
  
  global_delta <- median(vapply(
    pilot_results,
    function(x) {
      x$final_proposal$delta_final
    },
    numeric(1)
  ))
  
  global_initial_values <- list(
    theta_prior_z_offset = setNames(
      apply(
        pooled_theta_standardized,
        2L,
        median
      ),
      theta_names
    ),
    
    beta = setNames(
      apply(
        pooled_beta,
        2L,
        median
      ),
      colnames(pooled_beta)
    ),
    
    tau2_scale_ratio =
      median(pooled_tau2_ratio),
    
    aux_tau_scale_ratio =
      median(pooled_aux_ratio),
    
    eta =
      median(pooled_eta),
    
    delta_init = global_delta,
    kappa_cov_init = global_kappa,
    
    proposal_cor_init =
      pooled_theta_correlation,
    
    joint_covariance_standardized =
      global_joint_covariance_standardized,
    
    joint_theta_scale_init =
      JOINT_THETA_SCALE_INIT_PRODUCTION,
    
    pooled_tail_draws =
      nrow(pooled_theta_standardized)
  )
  
  parameter_vectors <- list(
    z_log_ell_cov =
      pooled_theta_standardized[
        ,
        "log_ell_cov"
      ],
    z_log_sigma_cov =
      pooled_theta_standardized[
        ,
        "log_sigma_cov"
      ],
    raw_log_ell_cov =
      pooled_theta_raw[
        ,
        "log_ell_cov"
      ],
    raw_log_sigma_cov =
      pooled_theta_raw[
        ,
        "log_sigma_cov"
      ],
    beta0 =
      pooled_beta[
        ,
        "beta0"
      ],
    tau2_cov_scale_ratio =
      pooled_tau2_ratio,
    aux_cov_scale_ratio =
      pooled_aux_ratio,
    eta =
      pooled_eta
  )
  
  gamma_columns <- setdiff(
    colnames(pooled_beta),
    "beta0"
  )
  
  for (column_name in gamma_columns) {
    parameter_vectors[[column_name]] <-
      pooled_beta[, column_name]
  }
  
  distribution_rows <- lapply(
    names(parameter_vectors),
    function(parameter_name) {
      values <-
        parameter_vectors[[parameter_name]]
      
      statistics <- summary_stats(values)
      
      data.frame(
        parameter = parameter_name,
        mean =
          unname(statistics["mean"]),
        sd =
          unname(statistics["sd"]),
        q025 =
          unname(statistics["q025"]),
        median =
          unname(statistics["median"]),
        q975 =
          unname(statistics["q975"]),
        pooled_tail_draws =
          length(values),
        stringsAsFactors = FALSE
      )
    }
  )
  
  per_pilot_rows <- lapply(
    seq_along(pilot_results),
    function(pilot_number) {
      x <- pilot_results[[pilot_number]]
      hp <- x$active_priors$hyperprior
      
      data.frame(
        pilot = pilot_number,
        covariate_seed =
          x$covariate_seed,
        mcmc_seed =
          x$mcmc_seed,
        
        prior_log_ell_mean =
          hp$log_ell_mean,
        prior_log_ell_sd =
          hp$log_ell_sd,
        prior_log_sigma_mean =
          hp$log_sig_mean,
        prior_log_sigma_sd =
          hp$log_sig_sd,
        tau_scale_cov =
          x$active_priors$tau_scale,
        
        log_ell_median = median(
          x$tail$theta_log[
            ,
            "log_ell_cov"
          ]
        ),
        log_sigma_median = median(
          x$tail$theta_log[
            ,
            "log_sigma_cov"
          ]
        ),
        beta0_median = median(
          x$tail$beta[
            ,
            "beta0"
          ]
        ),
        eta_median = median(
          x$tail$eta
        ),
        tau2_cov_ratio_median = median(
          x$tail$tau2_cov /
            x$active_priors$tau_scale^2
        ),
        aux_cov_ratio_median = median(
          x$tail$aux_tau_cov *
            x$active_priors$tau_scale^2
        ),
        
        delta_final =
          x$final_proposal$delta_final,
        kappa_cov_final =
          x$final_proposal$kappa_cov_final,
        rho_cov_final =
          x$final_proposal$proposal_cor_final[
            1L,
            2L
          ],
        joint_theta_scale_final =
          x$final_proposal$joint_theta_scale_final,
        
        blocked_acceptance_f =
          x$acceptance$blocked_f,
        blocked_acceptance_theta =
          x$acceptance$blocked_theta,
        joint_tune_acceptance =
          x$acceptance$joint_tune,
        runtime_sec =
          x$runtime_sec,
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  pilot_prior_summary <- do.call(
    rbind,
    lapply(
      pilot_results,
      `[[`,
      "prior_calibration_summary"
    )
  )
  
  row.names(pilot_prior_summary) <- NULL
  
  list(
    settings_signature =
      PILOT_SETTINGS_SIGNATURE,
    created_at =
      timestamp_now(),
    global_initial_values =
      global_initial_values,
    distribution_summary =
      do.call(
        rbind,
        distribution_rows
      ),
    per_pilot_summary =
      do.call(
        rbind,
        per_pilot_rows
      ),
    pilot_prior_calibration_summary =
      pilot_prior_summary,
    pilot_seed_table =
      pilot_seed_table
  )
}

load_or_build_global_calibration <- function() {
  if (isTRUE(RESUME_PILOTS) &&
      !isTRUE(FORCE_FRESH_PILOTS) &&
      file.exists(GLOBAL_CALIBRATION_FILE)) {
    candidate <- tryCatch(
      readRDS(GLOBAL_CALIBRATION_FILE),
      error = function(e) NULL
    )
    
    if (is.list(candidate) &&
        identical(
          candidate$settings_signature,
          PILOT_SETTINGS_SIGNATURE
        )) {
      cat(
        "\nLoaded the completed five-pilot self-calibrated GP-Cov + PU run.\n"
      )
      
      return(candidate)
    }
  }
  
  pilot_results <- run_all_pilots()
  
  calibration <- build_global_calibration(
    pilot_results
  )
  
  saveRDS(
    calibration,
    GLOBAL_CALIBRATION_FILE
  )
  
  global_pilot_calibration <- calibration
  
  save(
    global_pilot_calibration,
    file = GLOBAL_CALIBRATION_RDATA
  )
  
  write.csv(
    calibration$distribution_summary,
    file.path(
      OUT_DIR,
      "GP_Cov_PU_global_pilot_distribution_summary.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    calibration$per_pilot_summary,
    file.path(
      OUT_DIR,
      "GP_Cov_PU_global_pilot_per_pilot_summary.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    calibration$pilot_prior_calibration_summary,
    file.path(
      OUT_DIR,
      "GP_Cov_PU_pilot_prior_calibration_summary.csv"
    ),
    row.names = FALSE
  )
  
  rm(pilot_results)
  invisible(gc())
  
  calibration
}

global_pilot_calibration <-
  load_or_build_global_calibration()

GLOBAL_CALIBRATION_HASH <-
  file_md5_or_na(
    GLOBAL_CALIBRATION_FILE
  )

write.csv(
  global_pilot_calibration$distribution_summary,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_global_pilot_distribution_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  global_pilot_calibration$per_pilot_summary,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_global_pilot_per_pilot_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  global_pilot_calibration$pilot_prior_calibration_summary,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_pilot_prior_calibration_summary.csv"
  ),
  row.names = FALSE
)

cat("\n================ GLOBAL PRIOR-RELATIVE INITIAL VALUES ================\n")
print(
  global_pilot_calibration$global_initial_values
)

## ============================================================================
## 11. GEOMETRY-SPECIFIC PRODUCTION INITIALIZATION AND EVALUATION
## ============================================================================

make_global_initialization <- function(
    bundle,
    calibration
) {
  prior <- bundle$prior_calibration
  
  theta_names <- c(
    "log_ell_cov",
    "log_sigma_cov"
  )
  
  theta_initial <- prior$theta_prior_mean[
    theta_names
  ] + prior$theta_prior_sd[
    theta_names
  ] * calibration$theta_prior_z_offset[
    theta_names
  ]
  
  beta_initial <- as.numeric(
    calibration$beta
  )
  
  tau2_initial <- prior$tau_scale^2 *
    calibration$tau2_scale_ratio
  
  aux_initial <- calibration$aux_tau_scale_ratio /
    prior$tau_scale^2
  
  eta_initial <- as.numeric(
    calibration$eta
  )
  
  if (length(theta_initial) != 2L ||
      any(!is.finite(theta_initial)) ||
      length(beta_initial) != 11L ||
      any(!is.finite(beta_initial)) ||
      !is.finite(tau2_initial) ||
      tau2_initial <= 0 ||
      !is.finite(aux_initial) ||
      aux_initial <= 0 ||
      length(eta_initial) != 1L ||
      !is.finite(eta_initial) ||
      eta_initial <= 0 ||
      eta_initial >= 1) {
    stop(
      "The geometry-specific global initialization is invalid."
    )
  }
  
  names(beta_initial) <- names(
    calibration$beta
  )
  
  D_prior <- diag(
    prior$theta_prior_sd[theta_names],
    2L
  )
  
  joint_covariance_current <-
    sanitize_covariance_matrix(
      D_prior %*%
        calibration$joint_covariance_standardized %*%
        D_prior
    )
  
  dimnames(joint_covariance_current) <- list(
    theta_names,
    theta_names
  )
  
  list(
    state = list(
      theta_log = setNames(
        as.numeric(theta_initial),
        theta_names
      ),
      beta = beta_initial,
      tau2_cov =
        as.numeric(tau2_initial),
      aux_tau_cov =
        as.numeric(aux_initial),
      eta =
        as.numeric(eta_initial),
      ## Never transfer f between data sets.
      f = NULL
    ),
    
    joint_cov_base_init =
      joint_covariance_current
  )
}


evaluate_gp_cov_pu <- function(
    fit,
    replicate_number,
    runtime_sec,
    observed_y,
    hidden_positive_idx,
    signal_type
) {
  probabilities <- validate_probability_vector(
    fit$posterior_mean_scores,
    N_NODES,
    METHOD_NAME
  )
  
  predictive_T_draws <- validate_predictive_label_draws(
    predictive_T_draws = fit$predictive_T_draws,
    n_draws_expected = fit$n_saved,
    n_nodes_expected = N_NODES,
    method_label = METHOD_NAME
  )
  
  validate_sar_bundle_labels(
    observed_y = observed_y,
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = which(observed_y == 1L),
    signal_type = signal_type,
    context_label = paste0(
      "GP-Cov + PU SAR production replication ",
      replicate_number
    )
  )

  observed_y <- as.integer(observed_y)
  hidden_positive_idx <- sort(
    as.integer(hidden_positive_idx)
  )

  unlabeled_idx <- which(observed_y == 0L)
  hidden_truth_unlabeled <- as.integer(TRUE_T[unlabeled_idx] == 1L)
  
  overall_rank <- rank(
    -probabilities,
    ties.method = "average"
  )
  
  unlabeled_probabilities <- probabilities[unlabeled_idx]
  
  unlabeled_rank <- rank(
    -unlabeled_probabilities,
    ties.method = "average"
  )
  
  hidden_at_15 <- sum(
    hidden_truth_unlabeled == 1L &
      unlabeled_rank <= HIDDEN_TOP_K
  )
  
  overall_true_in_top_45 <- sum(
    TRUE_T == 1L &
      overall_rank <= OVERALL_RECALL_K
  )
  
  predictive_80 <- predictive_label_metrics(
    predictive_T_draws = predictive_T_draws,
    truth = TRUE_T,
    hidden_idx = hidden_positive_idx,
    level = 0.80
  )
  
  predictive_95 <- predictive_label_metrics(
    predictive_T_draws = predictive_T_draws,
    truth = TRUE_T,
    hidden_idx = hidden_positive_idx,
    level = 0.95
  )
  
  eta_q10 <- unname(
    quantile(fit$eta_draws, 0.10)
  )
  eta_q50 <- unname(
    quantile(fit$eta_draws, 0.50)
  )
  eta_q90 <- unname(
    quantile(fit$eta_draws, 0.90)
  )
  eta_q025 <- unname(
    quantile(fit$eta_draws, 0.025)
  )
  eta_q975 <- unname(
    quantile(fit$eta_draws, 0.975)
  )
  
  metrics <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    model_key = MODEL_KEY,
    method = METHOD_NAME,
    
    Overall_AUC = roc_auc(
      probabilities,
      TRUE_T
    ),
    
    Hidden_AUC = roc_auc(
      unlabeled_probabilities,
      hidden_truth_unlabeled
    ),
    
    Hidden_at_15 = hidden_at_15,
    
    Hidden_Recall_at_15 = hidden_at_15 / N_HIDDEN,
    
    Overall_Recall_at_45 = overall_true_in_top_45 / N_TRUE_POS,
    
    Overall_LogLoss = overall_logloss(
      probabilities,
      TRUE_T
    ),
    
    Hidden_LogLoss = hidden_logloss(
      probabilities,
      hidden_positive_idx
    ),
    
    Overall_Brier = overall_brier(
      probabilities,
      TRUE_T
    ),
    
    Hidden_Brier = hidden_brier(
      probabilities,
      hidden_positive_idx
    ),
    
    Overall_Coverage_80 = predictive_80$overall_coverage,
    Hidden_Coverage_80 = predictive_80$hidden_coverage,
    Overall_Interval_Score_80 =
      predictive_80$overall_interval_score,
    Hidden_Interval_Score_80 =
      predictive_80$hidden_interval_score,
    
    Overall_Coverage_95 = predictive_95$overall_coverage,
    Hidden_Coverage_95 = predictive_95$hidden_coverage,
    Overall_Interval_Score_95 =
      predictive_95$overall_interval_score,
    Hidden_Interval_Score_95 =
      predictive_95$hidden_interval_score,
    
    eta_posterior_mean = mean(fit$eta_draws),
    eta_posterior_sd = sd(fit$eta_draws),
    eta_q10 = eta_q10,
    eta_q50 = eta_q50,
    eta_q90 = eta_q90,
    eta_q025 = eta_q025,
    eta_q975 = eta_q975,
    
    anchor_max_mean = mean(fit$max_probability_draws),
    anchor_max_median = median(fit$max_probability_draws),
    anchor_max_q025 = unname(
      quantile(fit$max_probability_draws, 0.025)
    ),
    anchor_max_q975 = unname(
      quantile(fit$max_probability_draws, 0.975)
    ),
    anchor_Pr_over_090 = mean(
      fit$max_probability_draws > 0.90
    ),
    anchor_Pr_over_095 = mean(
      fit$max_probability_draws > 0.95
    ),
    anchor_Pr_over_099 = mean(
      fit$max_probability_draws > 0.99
    ),
    anchor_q99_median = median(
      fit$q99_probability_draws
    ),
    anchor_top5_median = median(
      fit$top5_mean_probability_draws
    ),
    beta0_posterior_mean = mean(fit$beta0_draws),
    tau_cov_posterior_mean = mean(fit$tau_cov_draws),
    ell_cov_posterior_mean = mean(
      fit$theta_draws[, "ell_cov"]
    ),
    sigma_cov_posterior_mean = mean(
      fit$theta_draws[, "sigma_cov"]
    ),
    acc_blocked_f = fit$acceptance$blocked_f,
    acc_blocked_theta = fit$acceptance$blocked_theta,
    acc_joint_tune = fit$acceptance$joint_tune,
    acc_joint_sampling = fit$acceptance$joint_sampling,
    acc_f_refresh_sampling = fit$acceptance$f_refresh_sampling,
    eta_moved_sampling = fit$acceptance$eta_moved_sampling,
    delta_final = fit$proposal$delta_final,
    kappa_cov_final = fit$proposal$kappa_cov_final,
    joint_theta_scale_final = fit$proposal$joint_theta_scale_final,
    runtime_sec = runtime_sec,
    stringsAsFactors = FALSE
  )
  
  node_scores <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    model_key = MODEL_KEY,
    method = METHOD_NAME,
    node = seq_len(N_NODES),
    T = TRUE_T,
    Y = observed_y,
    is_hidden_positive = as.integer(
      seq_len(N_NODES) %in% hidden_positive_idx
    ),
    positive_signal_group =
      POSITIVE_SIGNAL_GROUP,
    covariate_signal = as.integer(
      seq_len(N_NODES) %in% COVARIATE_SIGNAL_IDX
    ),
    network_signal = as.integer(
      seq_len(N_NODES) %in% NETWORK_SIGNAL_IDX
    ),
    signal_type = signal_type,
    probability = probabilities,
    overall_rank = overall_rank,
    unlabeled_rank = NA_real_,
    
    PI80_lower = predictive_80$lower,
    PI80_upper = predictive_80$upper,
    covered_80 = predictive_80$covered,
    interval_score_80 = predictive_80$score,
    
    PI95_lower = predictive_95$lower,
    PI95_upper = predictive_95$upper,
    covered_95 = predictive_95$covered,
    interval_score_95 = predictive_95$score,
    
    stringsAsFactors = FALSE
  )
  
  node_scores$unlabeled_rank[unlabeled_idx] <- unlabeled_rank
  
  list(
    metrics = metrics,
    node_scores = node_scores
  )
}

## ============================================================================
## 12. PRODUCTION SETTINGS SIGNATURE
## ============================================================================

PRODUCTION_SETTINGS_SIGNATURE <- list(
  script_version = SCRIPT_VERSION,
  evaluation_version = EVALUATION_VERSION,
  master_seed = MASTER_SEED,
  n_rep = N_REP,
  production_burnin =
    PRODUCTION_BURNIN,
  production_sampling =
    PRODUCTION_SAMPLING,
  production_thin =
    PRODUCTION_THIN,
  production_blocked_frac =
    PRODUCTION_BLOCKED_FRAC,
  global_calibration_hash =
    GLOBAL_CALIBRATION_HASH,
  use_saved_competing_geometries =
    USE_SAVED_COMPETING_GEOMETRIES,

  data_design = list(
    true_t = TRUE_T,
    true_zero_idx = TRUE_ZERO_IDX,
    true_positive_idx = TRUE_POS_IDX,
    covariates_only_idx = COV_ONLY_IDX,
    network_only_idx = NET_ONLY_IDX,
    both_idx = BOTH_IDX,
    covariate_signal_idx =
      COVARIATE_SIGNAL_IDX,
    network_signal_idx =
      NETWORK_SIGNAL_IDX,
    n_hidden_per_replication =
      N_HIDDEN,
    n_observed_positive_per_replication =
      N_OBSERVED_POS
  ),

  label_mechanism =
    PILOT_SETTINGS_SIGNATURE$label_mechanism,

  covariate_dgp =
    PILOT_SETTINGS_SIGNATURE$covariate_dgp,

  prior_calibration_rules =
    PRIOR_CALIBRATION_RULES,

  sampler_controls = list(
    adapt_block_cor =
      ADAPT_BLOCK_CORRELATION_PRODUCTION,
    targets = c(
      f = TARGET_F,
      theta = TARGET_THETA,
      joint = TARGET_JOINT
    ),
    max_log_sigma_cov =
      MAX_LOG_SIGMA_COV,
    min_log_sigma =
      MIN_LOG_SIGMA,
    logl_cov_bounds =
      LOGL_COV_BOUNDS,
    eta_slice_width =
      ETA_SLICE_WIDTH
  ),

  evaluation = list(
    version = EVALUATION_VERSION,
    hidden_top_k = HIDDEN_TOP_K,
    overall_recall_k =
      OVERALL_RECALL_K,
    predictive_interval_levels =
      PREDICTIVE_INTERVAL_LEVELS,
    predictive_label_distribution = paste0(
      "For every unit and posterior draw, use ",
      "T_i ~ Bernoulli(plogis(f_i))."
    ),
    interval_quantiles =
      "central discrete empirical quantiles with type=1",
    ground_truth =
      "fixed TRUE_T; hidden and observed positives both equal 1",
    H_intervals = FALSE,
    poisson_binomial = FALSE
  )
)

## ============================================================================
## 13. ONE PRODUCTION REPLICATION
## ============================================================================

run_one_production_replication <- function(
    replicate_number
) {
  checkpoint_file <- file.path(
    PRODUCTION_CHECKPOINT_DIR,
    sprintf(
      "GP_Cov_PU_replicate_%03d.rds",
      replicate_number
    )
  )
  
  if (isTRUE(RESUME_PRODUCTION) &&
      !isTRUE(FORCE_FRESH_PRODUCTION) &&
      file.exists(checkpoint_file)) {
    candidate <- tryCatch(
      readRDS(checkpoint_file),
      error = function(e) NULL
    )
    
    if (is.list(candidate) &&
        identical(
          candidate$settings_signature,
          PRODUCTION_SETTINGS_SIGNATURE
        ) &&
        identical(
          as.integer(candidate$replicate),
          as.integer(replicate_number)
        ) &&
        isTRUE(candidate$complete)) {
      cat(sprintf(
        paste0(
          "[%s] Reusing self-calibrated ",
          "GP-Cov + PU SAR replication %03d / %03d\n"
        ),
        timestamp_now(),
        replicate_number,
        N_REP
      ))
      
      return(candidate)
    }
  }
  
  replication_start <- Sys.time()
  
  bundle <- load_or_generate_production_bundle(
    replicate_number
  )
  
  observed_y <- as.integer(
    bundle$observed_y
  )

  hidden_positive_idx <- sort(
    as.integer(
      bundle$hidden_positive_idx
    )
  )

  observed_positive_idx <- sort(
    as.integer(
      bundle$observed_positive_idx
    )
  )

  signal_type <- as.character(
    bundle$signal_type
  )

  validate_sar_bundle_labels(
    observed_y = observed_y,
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = observed_positive_idx,
    signal_type = signal_type,
    context_label = paste0(
      "GP-Cov + PU SAR production replication ",
      replicate_number
    )
  )

  sar_selection_summary_now <-
    bundle$sar_selection_summary

  sar_selection_summary_now$replicate <-
    replicate_number

  sar_selection_summary_now <-
    sar_selection_summary_now[
      ,
      c(
        "scenario",
        "replicate",
        setdiff(
          names(sar_selection_summary_now),
          c("scenario", "replicate")
        )
      ),
      drop = FALSE
    ]

  sar_selection_diagnostics_now <-
    bundle$sar_selection_diagnostics

  sar_selection_diagnostics_now$replicate <-
    replicate_number

  sar_selection_diagnostics_now <-
    sar_selection_diagnostics_now[
      ,
      c(
        "scenario",
        "replicate",
        setdiff(
          names(sar_selection_diagnostics_now),
          c("scenario", "replicate")
        )
      ),
      drop = FALSE
    ]
  
  global_values <-
    global_pilot_calibration$global_initial_values
  
  initialization <- make_global_initialization(
    bundle,
    global_values
  )
  
  init_state <- initialization$state
  
  mcmc_seed <- as.integer(
    PRODUCTION_MCMC_SEED_BASE +
      10000L * replicate_number
  )
  
  predictive_seed <- as.integer(
    PRODUCTION_PREDICTIVE_SEED_BASE +
      10000L * replicate_number
  )
  
  cat("\n\n####################################################################\n")
  cat(sprintf(
    "# SELF-CALIBRATED GP-COV + PU SAR REPLICATION %d / %d\n",
    replicate_number,
    N_REP
  ))
  cat("####################################################################\n")
  cat("Geometry source       :",
      bundle$source, "\n")
  cat("Covariate seed        :",
      bundle$seeds$covariate_seed, "\n")
  cat("SAR seed              :",
      bundle$seeds$sar_seed, "\n")
  cat(
    "SAR hidden composition: cov=",
    sar_selection_summary_now$hidden_covariates_only,
    ", net=",
    sar_selection_summary_now$hidden_network_only,
    ", both=",
    sar_selection_summary_now$hidden_both,
    "\n",
    sep = ""
  )
  cat("Prior calibration time:",
      format_duration(
        bundle$prior_calibration$runtime_sec
      ),
      "\n")
  cat("MCMC seed             :",
      mcmc_seed, "\n")
  cat("Predictive-label seed :",
      predictive_seed, "\n")
  cat("Active kernel prior   :\n")
  print(
    bundle$prior_calibration$active_hyperprior
  )
  cat("Active tau scale      :",
      bundle$prior_calibration$tau_scale,
      "\n")
  
  fit_error <- NA_character_
  
  fit <- tryCatch({
    set.seed(mcmc_seed)
    
    run_gp_cov_pu_exact_chain(
      Y = observed_y,
      X_cov = bundle$X_cov,
      
      hyperprior_cov =
        bundle$prior_calibration$active_hyperprior,
      
      m_beta =
        beta0_prior_for_sampler$mean,
      s_beta =
        beta0_prior_for_sampler$sd,
      
      tau2_default =
        bundle$prior_calibration$tau2_default,
      scale_tau =
        bundle$prior_calibration$tau_scale,
      a_eta =
        eta_prior_for_sampler$a_eta,
      b_eta =
        eta_prior_for_sampler$b_eta,
      
      n_burnin =
        PRODUCTION_BURNIN,
      n_sample =
        PRODUCTION_SAMPLING,
      thin =
        PRODUCTION_THIN,
      blocked_frac =
        PRODUCTION_BLOCKED_FRAC,
      min_blocked_warmup =
        MIN_BLOCKED_WARMUP,
      
      init_state =
        init_state,
      
      delta_init =
        global_values$delta_init,
      kappa_cov_init =
        global_values$kappa_cov_init,
      proposal_cor_init =
        global_values$proposal_cor_init,
      
      adapt_block_cor =
        ADAPT_BLOCK_CORRELATION_PRODUCTION,
      cor_adapt_start =
        COR_ADAPT_START,
      cor_adapt_interval =
        COR_ADAPT_INTERVAL,
      cor_adapt_window =
        COR_ADAPT_WINDOW,
      cor_shrinkage =
        COR_SHRINKAGE,
      cor_max_abs =
        COR_MAX_ABS,
      
      joint_cov_base_init =
        initialization$joint_cov_base_init,
      joint_theta_scale_init =
        global_values$joint_theta_scale_init,
      f_refresh_after_joint =
        F_REFRESH_AFTER_JOINT,
      eta_slice_width =
        ETA_SLICE_WIDTH,
      initial_eta =
        INITIAL_ETA,
      
      nu_tau =
        NU_TAU,
      sample_tau =
        SAMPLE_TAU,
      
      adapt =
        ADAPT,
      adapt_start =
        ADAPT_START,
      adapt_interval =
        ADAPT_INTERVAL,
      target_f =
        TARGET_F,
      target_theta =
        TARGET_THETA,
      target_joint =
        TARGET_JOINT,
      adapt_rate_delta =
        ADAPT_RATE_DELTA,
      adapt_rate_theta =
        ADAPT_RATE_THETA,
      adapt_rate_joint_theta =
        ADAPT_RATE_JOINT_THETA,
      t0_adapt =
        T0_ADAPT,
      power_adapt =
        POWER_ADAPT,
      
      delta_min =
        DELTA_MIN,
      delta_max =
        DELTA_MAX,
      kappa_min =
        KAPPA_MIN,
      kappa_max =
        KAPPA_MAX,
      joint_theta_scale_min =
        JOINT_THETA_SCALE_MIN,
      joint_theta_scale_max =
        JOINT_THETA_SCALE_MAX,
      
      emp_cov_tail =
        EMP_COV_TAIL,
      cov_jitter =
        COV_JITTER,
      max_log_sigma_cov =
        MAX_LOG_SIGMA_COV,
      min_log_sigma =
        MIN_LOG_SIGMA,
      logl_cov_bounds =
        LOGL_COV_BOUNDS,
      jitter =
        JITTER,
      
      state_tail_length =
        0L,
      predictive_seed =
        predictive_seed,
      
      phase_label = paste0(
        "GP-Cov + PU SAR replication ",
        sprintf(
          "%03d",
          replicate_number
        )
      ),
      burnin_progress_every =
        PRODUCTION_BURNIN_PROGRESS_EVERY,
      sample_progress_every =
        PRODUCTION_SAMPLE_PROGRESS_EVERY,
      verbose = TRUE
    )
  }, error = function(e) {
    fit_error <<- conditionMessage(e)
    NULL
  })
  
  runtime_sec <- as.numeric(difftime(
    Sys.time(),
    replication_start,
    units = "secs"
  ))
  
  if (is.null(fit)) {
    failure <- data.frame(
      scenario = SCENARIO_NAME,
      replicate = replicate_number,
      model_key = MODEL_KEY,
      method = METHOD_NAME,
      runtime_sec = runtime_sec,
      error_message = fit_error,
      stringsAsFactors = FALSE
    )
    
    result <- list(
      settings_signature =
        PRODUCTION_SETTINGS_SIGNATURE,
      replicate = replicate_number,
      complete = FALSE,
      failure = failure,
      prior_calibration_summary =
        bundle$prior_calibration$summary,
      sar_selection_summary =
        sar_selection_summary_now,
      sar_selection_diagnostics =
        sar_selection_diagnostics_now,
      label_design = list(
        observed_y = observed_y,
        hidden_positive_idx =
          hidden_positive_idx,
        observed_positive_idx =
          observed_positive_idx,
        signal_type = signal_type
      ),
      time_start = replication_start,
      time_end = Sys.time(),
      runtime_sec = runtime_sec
    )
    
    saveRDS(
      result,
      checkpoint_file
    )
    
    if (isTRUE(FAIL_IF_PRODUCTION_FAILS)) {
      stop(
        "GP-Cov + PU failed in replication ",
        replicate_number,
        ": ",
        fit_error,
        ". Partial checkpoint saved to ",
        checkpoint_file,
        "."
      )
    }
    
    return(result)
  }
  
  evaluation <- evaluate_gp_cov_pu(
    fit = fit,
    replicate_number =
      replicate_number,
    runtime_sec =
      runtime_sec,
    observed_y =
      observed_y,
    hidden_positive_idx =
      hidden_positive_idx,
    signal_type =
      signal_type
  )
  
  prior_summary <-
    bundle$prior_calibration$summary
  
  prior_summary$scenario <-
    SCENARIO_NAME
  prior_summary$replicate <-
    replicate_number
  prior_summary$geometry_source <-
    bundle$source
  
  prior_summary <- prior_summary[
    ,
    c(
      "scenario",
      "replicate",
      "geometry_source",
      setdiff(
        names(prior_summary),
        c(
          "scenario",
          "replicate",
          "geometry_source"
        )
      )
    ),
    drop = FALSE
  ]
  
  write.csv(
    prior_summary,
    file.path(
      PRIOR_CALIBRATION_DIR,
      sprintf(
        "GP_Cov_PU_replicate_%03d_prior_calibration.csv",
        replicate_number
      )
    ),
    row.names = FALSE
  )
  
  geometry_summary <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    geometry_source = bundle$source,
    covariate_seed =
      bundle$seeds$covariate_seed,
    sar_seed =
      bundle$seeds$sar_seed,
    lambda_SAR =
      LAMBDA_SAR,
    n = N_NODES,
    n_covariates =
      ncol(bundle$X_cov),
    n_true_zero =
      N_TRUE_ZERO,
    n_true_positive =
      N_TRUE_POS,
    n_observed_positive =
      length(observed_positive_idx),
    n_hidden_positive =
      length(hidden_positive_idx),
    hidden_covariates_only =
      sar_selection_summary_now$hidden_covariates_only,
    hidden_network_only =
      sar_selection_summary_now$hidden_network_only,
    hidden_both =
      sar_selection_summary_now$hidden_both,
    mean_oracle_score_hidden =
      sar_selection_summary_now$mean_score_hidden,
    mean_oracle_score_observed =
      sar_selection_summary_now$mean_score_observed,
    observed_minus_hidden_score =
      sar_selection_summary_now$observed_minus_hidden_score,
    n_covariate_signal =
      length(COVARIATE_SIGNAL_IDX),
    n_network_signal =
      length(NETWORK_SIGNAL_IDX),
    n_both_signal =
      length(intersect(
        COVARIATE_SIGNAL_IDX,
        NETWORK_SIGNAL_IDX
      )),
    prior_calibration_runtime_sec =
      bundle$prior_calibration$runtime_sec,
    tau_scale_cov =
      bundle$prior_calibration$tau_scale,
    stringsAsFactors = FALSE
  )
  
  result <- list(
    settings_signature =
      PRODUCTION_SETTINGS_SIGNATURE,
    replicate = replicate_number,
    complete = TRUE,
    
    metrics =
      evaluation$metrics,
    node_scores =
      evaluation$node_scores,
    
    prior_calibration_summary =
      prior_summary,
    
    active_priors = list(
      hyperprior =
        bundle$prior_calibration$active_hyperprior,
      hyperprior_unrounded =
        bundle$prior_calibration$hyperprior_unrounded,
      tau_scale =
        bundle$prior_calibration$tau_scale,
      tau_scale_unrounded =
        bundle$prior_calibration$tau_scale_unrounded,
      tau2_default =
        bundle$prior_calibration$tau2_default,
      theta_prior_mean =
        bundle$prior_calibration$theta_prior_mean,
      theta_prior_sd =
        bundle$prior_calibration$theta_prior_sd
    ),
    
    geometry_summary =
      geometry_summary,

    sar_selection_summary =
      sar_selection_summary_now,
    sar_selection_diagnostics =
      sar_selection_diagnostics_now,
    label_design = list(
      observed_y = observed_y,
      hidden_positive_idx =
        hidden_positive_idx,
      observed_positive_idx =
        observed_positive_idx,
      signal_type = signal_type
    ),

    final_state =
      fit$state,
    final_proposal =
      fit$proposal,
    acceptance =
      fit$acceptance,
    n_saved =
      fit$n_saved,
    
    mcmc_seed =
      mcmc_seed,
    predictive_label_seed =
      predictive_seed,
    
    time_start =
      replication_start,
    time_end =
      Sys.time(),
    runtime_sec =
      runtime_sec,
    
    failure = NULL
  )
  
  saveRDS(
    result,
    checkpoint_file
  )
  
  rm(
    bundle,
    initialization,
    init_state,
    fit,
    evaluation,
    observed_y,
    hidden_positive_idx,
    observed_positive_idx,
    signal_type
  )
  
  invisible(gc())
  
  result
}


## ============================================================================
## 14. RUN THE PRODUCTION REPLICATIONS
## ============================================================================

if (
  length(PRODUCTION_REPLICATES) > 0L &&
  any(!PRODUCTION_REPLICATES %in% seq_len(N_REP))
) {
  stop(
    "PRODUCTION_REPLICATES contains indices outside 1:N_REP."
  )
}

if (anyDuplicated(PRODUCTION_REPLICATES)) {
  stop(
    "PRODUCTION_REPLICATES contains duplicated replication indices."
  )
}

production_results <- vector(
  "list",
  N_REP
)

production_start <- Sys.time()
n_loaded_complete <- 0L
n_computed_now <- 0L

for (
  batch_position in seq_along(
    PRODUCTION_REPLICATES
  )
) {
  replicate_number <- as.integer(
    PRODUCTION_REPLICATES[batch_position]
  )
  
  checkpoint_file <- file.path(
    PRODUCTION_CHECKPOINT_DIR,
    sprintf(
      "GP_Cov_PU_replicate_%03d.rds",
      replicate_number
    )
  )
  
  was_complete_before <- FALSE
  
  if (
    isTRUE(RESUME_PRODUCTION) &&
    !isTRUE(FORCE_FRESH_PRODUCTION) &&
    file.exists(checkpoint_file)
  ) {
    prior_checkpoint <- tryCatch(
      readRDS(checkpoint_file),
      error = function(e) NULL
    )
    
    was_complete_before <-
      is.list(prior_checkpoint) &&
      identical(
        prior_checkpoint$settings_signature,
        PRODUCTION_SETTINGS_SIGNATURE
      ) &&
      identical(
        as.integer(prior_checkpoint$replicate),
        replicate_number
      ) &&
      isTRUE(prior_checkpoint$complete)
  }
  
  production_results[[replicate_number]] <-
    run_one_production_replication(
      replicate_number
    )
  
  if (
    !isTRUE(
      production_results[[replicate_number]]$complete
    )
  ) {
    stop(
      "Production replication ",
      replicate_number,
      " did not complete."
    )
  }
  
  if (was_complete_before) {
    n_loaded_complete <-
      n_loaded_complete + 1L
  } else {
    n_computed_now <-
      n_computed_now + 1L
  }
  
  elapsed_sec <- as.numeric(
    difftime(
      Sys.time(),
      production_start,
      units = "secs"
    )
  )
  
  projected_total_sec <-
    elapsed_sec /
    batch_position *
    length(PRODUCTION_REPLICATES)
  
  cat(sprintf(
    paste0(
      "[%s] Completed GP-Cov + PU replication ",
      "%03d/%03d | batch %02d/%02d | ",
      "elapsed %s | projected batch total %s\n"
    ),
    timestamp_now(),
    replicate_number,
    N_REP,
    batch_position,
    length(PRODUCTION_REPLICATES),
    format_duration(elapsed_sec),
    format_duration(projected_total_sec)
  ))
}

production_end <- Sys.time()

## Worker runs stop here after saving their own checkpoints.
if (isTRUE(BATCH_ONLY)) {
  batch_description <- if (
    length(PRODUCTION_REPLICATES) > 0L
  ) {
    paste0(
      min(PRODUCTION_REPLICATES),
      ":",
      max(PRODUCTION_REPLICATES)
    )
  } else {
    "none"
  }
  
  cat(
    "\nCompleted GP-Cov + PU production batch: ",
    batch_description,
    "\nAll replication checkpoints were saved in:\n",
    normalizePath(
      PRODUCTION_CHECKPOINT_DIR,
      winslash = "/",
      mustWork = FALSE
    ),
    "\n",
    sep = ""
  )
  
  quit(save = "no")
}

## Combine-only mode: load all 100 completed checkpoints.
for (
  replicate_number in seq_len(N_REP)
) {
  checkpoint_file <- file.path(
    PRODUCTION_CHECKPOINT_DIR,
    sprintf(
      "GP_Cov_PU_replicate_%03d.rds",
      replicate_number
    )
  )
  
  if (!file.exists(checkpoint_file)) {
    stop(
      "Missing checkpoint for replication ",
      replicate_number,
      ": ",
      checkpoint_file
    )
  }
  
  candidate <- tryCatch(
    readRDS(checkpoint_file),
    error = function(e) NULL
  )
  
  if (
    !is.list(candidate) ||
    !identical(
      candidate$settings_signature,
      PRODUCTION_SETTINGS_SIGNATURE
    ) ||
    !identical(
      as.integer(candidate$replicate),
      as.integer(replicate_number)
    ) ||
    !isTRUE(candidate$complete)
  ) {
    stop(
      "Checkpoint for replication ",
      replicate_number,
      " is incomplete or belongs to different settings."
    )
  }
  
  production_results[[replicate_number]] <-
    candidate
}

n_loaded_complete <- N_REP
n_computed_now <- 0L


## ============================================================================
## 15. COMBINE PRODUCTION RESULTS
## ============================================================================

if (any(!vapply(
  production_results,
  function(x) {
    is.list(x) &&
      isTRUE(x$complete)
  },
  logical(1)
))) {
  stop(
    "Not all 100 GP-Cov + PU production replications completed."
  )
}

metrics_by_replicate <- do.call(
  rbind,
  lapply(
    production_results,
    `[[`,
    "metrics"
  )
)

row.names(metrics_by_replicate) <- NULL

scores_all_nodes <- do.call(
  rbind,
  lapply(
    production_results,
    `[[`,
    "node_scores"
  )
)

row.names(scores_all_nodes) <- NULL

geometry_summary <- do.call(
  rbind,
  lapply(
    production_results,
    `[[`,
    "geometry_summary"
  )
)

row.names(geometry_summary) <- NULL

sar_selection_summary <- do.call(
  rbind,
  lapply(
    production_results,
    `[[`,
    "sar_selection_summary"
  )
)

row.names(sar_selection_summary) <- NULL

sar_selection_diagnostics <- do.call(
  rbind,
  lapply(
    production_results,
    `[[`,
    "sar_selection_diagnostics"
  )
)

row.names(sar_selection_diagnostics) <- NULL

sar_selection_metric_names <- c(
  "hidden_covariates_only",
  "hidden_network_only",
  "hidden_both",
  "mean_score_hidden",
  "mean_score_observed",
  "observed_minus_hidden_score"
)

sar_selection_mc_summary <- do.call(
  rbind,
  lapply(
    sar_selection_metric_names,
    function(metric_name) {
      statistics <- summary_stats(
        sar_selection_summary[[metric_name]]
      )

      data.frame(
        scenario = SCENARIO_NAME,
        quantity = metric_name,
        mean = unname(statistics["mean"]),
        sd = unname(statistics["sd"]),
        mc_se = unname(statistics["mc_se"]),
        q025 = unname(statistics["q025"]),
        median = unname(statistics["median"]),
        q975 = unname(statistics["q975"]),
        stringsAsFactors = FALSE
      )
    }
  )
)


prior_calibration_summary <- do.call(
  rbind,
  lapply(
    production_results,
    `[[`,
    "prior_calibration_summary"
  )
)

row.names(prior_calibration_summary) <- NULL

failure_parts <- Filter(
  Negate(is.null),
  lapply(
    production_results,
    `[[`,
    "failure"
  )
)

failure_log <- if (length(failure_parts) == 0L) {
  data.frame(
    scenario = character(0),
    replicate = integer(0),
    model_key = character(0),
    method = character(0),
    runtime_sec = numeric(0),
    error_message = character(0),
    stringsAsFactors = FALSE
  )
} else {
  do.call(
    rbind,
    failure_parts
  )
}

if (nrow(metrics_by_replicate) != N_REP) {
  stop(
    "The metric table does not contain exactly 100 rows."
  )
}

if (nrow(scores_all_nodes) !=
    N_REP * N_NODES) {
  stop(
    "The node-score table does not contain exactly 100 x 400 rows."
  )
}

if (nrow(prior_calibration_summary) !=
    N_REP) {
  stop(
    "The prior-calibration table does not contain exactly 100 rows."
  )
}

if (
  nrow(sar_selection_summary) != N_REP ||
  nrow(sar_selection_diagnostics) !=
    N_REP * N_TRUE_POS
) {
  stop(
    "The SAR-selection tables do not contain the expected number of rows."
  )
}

metrics_by_replicate <-
  metrics_by_replicate[
    order(
      metrics_by_replicate$replicate
    ),
    ,
    drop = FALSE
  ]

scores_all_nodes <-
  scores_all_nodes[
    order(
      scores_all_nodes$replicate,
      scores_all_nodes$node
    ),
    ,
    drop = FALSE
  ]

geometry_summary <-
  geometry_summary[
    order(
      geometry_summary$replicate
    ),
    ,
    drop = FALSE
  ]

sar_selection_summary <-
  sar_selection_summary[
    order(
      sar_selection_summary$replicate
    ),
    ,
    drop = FALSE
  ]

sar_selection_diagnostics <-
  sar_selection_diagnostics[
    order(
      sar_selection_diagnostics$replicate,
      sar_selection_diagnostics$node
    ),
    ,
    drop = FALSE
  ]

prior_calibration_summary <-
  prior_calibration_summary[
    order(
      prior_calibration_summary$replicate
    ),
    ,
    drop = FALSE
  ]

write.csv(
  metrics_by_replicate,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_metrics_by_replicate.csv"
  ),
  row.names = FALSE
)

write.csv(
  scores_all_nodes,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_all_node_probabilities.csv"
  ),
  row.names = FALSE
)

write.csv(
  geometry_summary,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_geometry_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  sar_selection_summary,
  file.path(
    OUT_DIR,
    "SAR_selection_summary_by_replicate.csv"
  ),
  row.names = FALSE
)

write.csv(
  sar_selection_diagnostics,
  file.path(
    OUT_DIR,
    "SAR_selection_diagnostics_all_positives.csv"
  ),
  row.names = FALSE
)

write.csv(
  sar_selection_mc_summary,
  file.path(
    OUT_DIR,
    "SAR_selection_Monte_Carlo_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  prior_calibration_summary,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_prior_calibration_by_replicate.csv"
  ),
  row.names = FALSE
)

write.csv(
  failure_log,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_failure_log.csv"
  ),
  row.names = FALSE
)

prior_summary_variables <- c(
  "n_mean_features",
  "Z_rms_factor",
  "log_ell_mean_unrounded",
  "log_ell_sd_unrounded",
  "log_sig_mean_unrounded",
  "log_sig_sd_unrounded",
  "log_ell_mean_active",
  "log_ell_sd_active",
  "log_sig_mean_active",
  "log_sig_sd_active",
  "ell_median",
  "sigma_raw_median",
  "tau_scale_unrounded",
  "tau_scale_active",
  "prior_calibration_runtime_sec"
)

prior_calibration_distribution_summary <-
  do.call(
    rbind,
    lapply(
      prior_summary_variables,
      function(variable_name) {
        statistics <- summary_stats(
          prior_calibration_summary[[variable_name]]
        )
        
        data.frame(
          expert = "cov",
          parameter = variable_name,
          n_rep =
            nrow(prior_calibration_summary),
          mean =
            unname(statistics["mean"]),
          sd =
            unname(statistics["sd"]),
          mc_se =
            unname(statistics["mc_se"]),
          q025 =
            unname(statistics["q025"]),
          median =
            unname(statistics["median"]),
          q975 =
            unname(statistics["q975"]),
          stringsAsFactors = FALSE
        )
      }
    )
  )

row.names(
  prior_calibration_distribution_summary
) <- NULL

write.csv(
  prior_calibration_distribution_summary,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_prior_calibration_distribution_summary.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 16. MONTE CARLO SUMMARIES
## ============================================================================

main_metric_names <- c(
  "Overall_AUC",
  "Hidden_AUC",
  "Overall_LogLoss",
  "Hidden_LogLoss",
  "Overall_Brier",
  "Hidden_Brier",
  "Overall_Recall_at_45",
  "Hidden_Recall_at_15",
  "Overall_Coverage_80",
  "Hidden_Coverage_80",
  "Overall_Interval_Score_80",
  "Hidden_Interval_Score_80",
  "Overall_Coverage_95",
  "Hidden_Coverage_95",
  "Overall_Interval_Score_95",
  "Hidden_Interval_Score_95"
)

metric_direction <- c(
  Overall_AUC = "higher",
  Hidden_AUC = "higher",
  Overall_LogLoss = "lower",
  Hidden_LogLoss = "lower",
  Overall_Brier = "lower",
  Hidden_Brier = "lower",
  Overall_Recall_at_45 = "higher",
  Hidden_Recall_at_15 = "higher",
  Overall_Coverage_80 =
    "higher_interpret_with_interval_score",
  Hidden_Coverage_80 =
    "higher_interpret_with_interval_score",
  Overall_Interval_Score_80 = "lower",
  Hidden_Interval_Score_80 = "lower",
  Overall_Coverage_95 =
    "higher_interpret_with_interval_score",
  Hidden_Coverage_95 =
    "higher_interpret_with_interval_score",
  Overall_Interval_Score_95 = "lower",
  Hidden_Interval_Score_95 = "lower"
)

metric_summary_long <- do.call(
  rbind,
  lapply(
    main_metric_names,
    function(metric_name) {
      values <-
        metrics_by_replicate[[metric_name]]
      
      statistics <- summary_stats(values)
      
      data.frame(
        scenario = SCENARIO_NAME,
        method = METHOD_NAME,
        metric = metric_name,
        direction =
          unname(metric_direction[metric_name]),
        n_rep =
          sum(is.finite(values)),
        mean =
          unname(statistics["mean"]),
        sd =
          unname(statistics["sd"]),
        mc_se =
          unname(statistics["mc_se"]),
        q025 =
          unname(statistics["q025"]),
        median =
          unname(statistics["median"]),
        q975 =
          unname(statistics["q975"]),
        stringsAsFactors = FALSE
      )
    }
  )
)

row.names(metric_summary_long) <- NULL

write.csv(
  metric_summary_long,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_metric_summary_long.csv"
  ),
  row.names = FALSE
)

main_table_means <- data.frame(
  scenario = SCENARIO_NAME,
  method = METHOD_NAME,
  stringsAsFactors = FALSE
)

for (metric_name in main_metric_names) {
  main_table_means[[metric_name]] <- mean(
    metrics_by_replicate[[metric_name]],
    na.rm = TRUE
  )
}

write.csv(
  main_table_means,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_FINAL_MAIN_TABLE_means.csv"
  ),
  row.names = FALSE
)

format_mean_sd <- function(
    mean_value,
    sd_value,
    digits = 3L
) {
  if (!is.finite(mean_value)) {
    return("")
  }
  
  if (!is.finite(sd_value)) {
    return(formatC(
      mean_value,
      format = "f",
      digits = digits
    ))
  }
  
  paste0(
    formatC(
      mean_value,
      format = "f",
      digits = digits
    ),
    " (",
    formatC(
      sd_value,
      format = "f",
      digits = digits
    ),
    ")"
  )
}

paper_column_labels <- c(
  Overall_AUC = "Overall AUC",
  Hidden_AUC = "Hidden AUC",
  Overall_LogLoss = "Overall LogLoss",
  Hidden_LogLoss = "Hidden LogLoss",
  Overall_Brier = "Overall Brier",
  Hidden_Brier = "Hidden Brier",
  Overall_Recall_at_45 =
    "Overall Recall@45",
  Hidden_Recall_at_15 =
    "Hidden Recall@15",
  Overall_Coverage_80 =
    "Overall Coverage 80%",
  Hidden_Coverage_80 =
    "Hidden Coverage 80%",
  Overall_Interval_Score_80 =
    "Overall Interval Score 80%",
  Hidden_Interval_Score_80 =
    "Hidden Interval Score 80%",
  Overall_Coverage_95 =
    "Overall Coverage 95%",
  Hidden_Coverage_95 =
    "Hidden Coverage 95%",
  Overall_Interval_Score_95 =
    "Overall Interval Score 95%",
  Hidden_Interval_Score_95 =
    "Hidden Interval Score 95%"
)

paper_table_mean_sd <- data.frame(
  Method = METHOD_NAME,
  stringsAsFactors = FALSE
)

for (metric_name in main_metric_names) {
  values <-
    metrics_by_replicate[[metric_name]]
  
  paper_table_mean_sd[[paper_column_labels[metric_name]]] <- format_mean_sd(
    mean(values, na.rm = TRUE),
    stats::sd(values, na.rm = TRUE),
    digits = 3L
  )
}

write.csv(
  paper_table_mean_sd,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_FINAL_PAPER_TABLE_mean_SD.csv"
  ),
  row.names = FALSE,
  na = ""
)

uncertainty_metric_names <- c(
  "Overall_Coverage_80",
  "Hidden_Coverage_80",
  "Overall_Interval_Score_80",
  "Hidden_Interval_Score_80",
  "Overall_Coverage_95",
  "Hidden_Coverage_95",
  "Overall_Interval_Score_95",
  "Hidden_Interval_Score_95"
)

uncertainty_table_mean_sd <-
  paper_table_mean_sd[
    ,
    c(
      "Method",
      unname(
        paper_column_labels[
          uncertainty_metric_names
        ]
      )
    ),
    drop = FALSE
  ]

write.csv(
  uncertainty_table_mean_sd,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_FINAL_UNCERTAINTY_TABLE_mean_SD.csv"
  ),
  row.names = FALSE,
  na = ""
)

additional_metric_names <- c(
  "Hidden_at_15",
  "eta_posterior_mean",
  "eta_posterior_sd",
  "eta_q10",
  "eta_q50",
  "eta_q90",
  "eta_q025",
  "eta_q975",
  "anchor_max_median",
  "anchor_Pr_over_090",
  "anchor_Pr_over_095",
  "anchor_Pr_over_099",
  "anchor_q99_median",
  "anchor_top5_median",
  "beta0_posterior_mean",
  "tau_cov_posterior_mean",
  "ell_cov_posterior_mean",
  "sigma_cov_posterior_mean",
  "acc_blocked_f",
  "acc_blocked_theta",
  "acc_joint_tune",
  "acc_joint_sampling",
  "acc_f_refresh_sampling",
  "eta_moved_sampling",
  "delta_final",
  "kappa_cov_final",
  "joint_theta_scale_final",
  "runtime_sec"
)

additional_diagnostic_summary <- do.call(
  rbind,
  lapply(
    additional_metric_names,
    function(metric_name) {
      values <-
        metrics_by_replicate[[metric_name]]
      
      statistics <- summary_stats(values)
      
      data.frame(
        method = METHOD_NAME,
        metric = metric_name,
        n_rep =
          sum(is.finite(values)),
        mean =
          unname(statistics["mean"]),
        sd =
          unname(statistics["sd"]),
        mc_se =
          unname(statistics["mc_se"]),
        q025 =
          unname(statistics["q025"]),
        median =
          unname(statistics["median"]),
        q975 =
          unname(statistics["q975"]),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  additional_diagnostic_summary,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_additional_diagnostic_summary.csv"
  ),
  row.names = FALSE
)

eta_diagnostics_by_replicate <- metrics_by_replicate[
  ,
  c(
    "replicate",
    "eta_posterior_mean",
    "eta_posterior_sd",
    "eta_q10",
    "eta_q50",
    "eta_q90",
    "eta_q025",
    "eta_q975"
  ),
  drop = FALSE
]

write.csv(
  eta_diagnostics_by_replicate,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_eta_descriptive_diagnostics_by_replicate.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 17. SAVE THE COMPLETE OUTPUT BUNDLE
## ============================================================================

total_production_runtime_sec <- sum(
  vapply(
    production_results,
    function(x) {
      as.numeric(x$runtime_sec)
    },
    numeric(1)
  ),
  na.rm = TRUE
)

run_config <- list(
  script_version =
    SCRIPT_VERSION,
  evaluation_version =
    EVALUATION_VERSION,
  created_at =
    timestamp_now(),
  scenario =
    SCENARIO_NAME,
  method =
    METHOD_NAME,
  
  model_definition = paste0(
    "GP-Cov + PU: f = beta0 + X gamma_cov + g_cov, ",
    "g_cov ~ GP(0,K_RBF), ",
    "Y|f,eta ~ Bernoulli((1-eta)*plogis(f))."
  ),
  
  prior_calibration_mode =
    "exact_per_covariate_geometry",
  prior_calibration_rules =
    PRIOR_CALIBRATION_RULES,
  beta0_prior =
    beta0_prior_for_sampler,
  eta_prior =
    eta_prior_for_sampler,
  fixed_tau_cov_scale =
    TAU_SCALE_COV_ACTIVE,
  fixed_tau_cov_scale_unrounded =
    TAU_SCALE_COV_UNROUNDED,
  active_prior_rounding =
    ROUND_ACTIVE_PRIORS,
  active_prior_round_digits =
    PRIOR_ROUND_DIGITS,
  
  global_calibration_file =
    GLOBAL_CALIBRATION_FILE,
  global_calibration_hash =
    GLOBAL_CALIBRATION_HASH,
  
  n_pilot =
    N_PILOT,
  pilot_burnin =
    PILOT_BURNIN,
  pilot_tail_length =
    PILOT_TAIL_LENGTH,
  
  n_rep =
    N_REP,
  production_burnin =
    PRODUCTION_BURNIN,
  production_sampling =
    PRODUCTION_SAMPLING,
  production_thin =
    PRODUCTION_THIN,
  
  true_labels =
    TRUE_T,
  positive_signal_group =
    POSITIVE_SIGNAL_GROUP,

  label_mechanism = list(
    type = "covariate-driven SAR",
    lambda = LAMBDA_SAR,
    z_clip = SAR_Z_CLIP,
    n_hidden_per_replication =
      N_HIDDEN,
    n_observed_positive_per_replication =
      N_OBSERVED_POS,
    selection = paste0(
      "Exactly 15 of the 45 true positives are sampled without replacement ",
      "with weights proportional to exp(-lambda * z_i), where z_i is the ",
      "oracle covariate fraud-likeness score standardized among the true ",
      "positives and clipped to [-3,3]."
    ),
    oracle_score_used_for_fitting =
      FALSE
  ),

  sar_selection_summary =
    sar_selection_summary,
  sar_selection_Monte_Carlo_summary =
    sar_selection_mc_summary,

  production_seed_table =
    production_seed_table,
  pilot_seed_table =
    pilot_seed_table,
  
  use_saved_competing_geometries =
    USE_SAVED_COMPETING_GEOMETRIES,
  
  n_loaded_complete =
    n_loaded_complete,
  n_computed_now =
    n_computed_now,
  total_production_runtime_sec =
    total_production_runtime_sec,
  
  predictive_interval_levels =
    PREDICTIVE_INTERVAL_LEVELS,
  predictive_label_construction = paste0(
    "For every unit and posterior draw, sample T_i from ",
    "Bernoulli(plogis(f_i)). Use central discrete empirical ",
    "quantiles with type=1."
  ),
  evaluation_truth = paste0(
    "Intervals are compared with fixed TRUE_T: 0 for true zeros and 1 for ",
    "both observed and hidden positives. The ground truth is never regenerated."
  ),
  H_intervals_used =
    FALSE,
  poisson_binomial_used =
    FALSE,
  application_dependency =
    "none"
)

saveRDS(
  run_config,
  file.path(
    OUT_DIR,
    "GP_Cov_PU_run_config.rds"
  )
)

save(
  global_pilot_calibration,
  production_results,
  metrics_by_replicate,
  scores_all_nodes,
  geometry_summary,
  sar_selection_summary,
  sar_selection_diagnostics,
  sar_selection_mc_summary,
  prior_calibration_summary,
  prior_calibration_distribution_summary,
  failure_log,
  metric_summary_long,
  main_table_means,
  paper_table_mean_sd,
  uncertainty_table_mean_sd,
  additional_diagnostic_summary,
  eta_diagnostics_by_replicate,
  production_seed_table,
  pilot_seed_table,
  TRUE_T,
  TRUE_ZERO_IDX,
  TRUE_POS_IDX,
  COV_ONLY_IDX,
  NET_ONLY_IDX,
  BOTH_IDX,
  COVARIATE_SIGNAL_IDX,
  NETWORK_SIGNAL_IDX,
  POSITIVE_SIGNAL_GROUP,
  LAMBDA_SAR,
  SAR_Z_CLIP,
  beta0_prior_for_sampler,
  eta_prior_for_sampler,
  PRIOR_CALIBRATION_RULES,
  TAU_SCALE_COV_ACTIVE,
  TAU_SCALE_COV_UNROUNDED,
  run_config,
  PILOT_SETTINGS_SIGNATURE,
  PRODUCTION_SETTINGS_SIGNATURE,
  file = file.path(
    OUT_DIR,
    "SIMULATION_100_GP_COV_PU_SELF_CALIBRATING_SAR_COMPLETE_RESULTS.RData"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    OUT_DIR,
    "sessionInfo_GP_Cov_PU_self_calibrating_SAR_simulation.txt"
  )
)

## ============================================================================
## 18. FINAL CONSOLE OUTPUT
## ============================================================================

cat("\n================ RUN ACCOUNTING ================\n")
cat("Pilot datasets          :", N_PILOT, "\n")
cat("Production replications :", N_REP, "\n")
cat("Computed in this run    :", n_computed_now, "\n")
cat("Loaded from checkpoint  :", n_loaded_complete, "\n")
cat("Production runtime      :",
    format_duration(total_production_runtime_sec), "\n")
cat("Prior calibration mode  : exact per pilot/replication\n")
cat("Application dependency  : none\n")
cat("SAR lambda              :", LAMBDA_SAR, "\n")

cat("\n================ SAR SELECTION MONTE CARLO SUMMARY ================\n")
print(
  sar_selection_mc_summary,
  row.names = FALSE
)

cat("\n================ FINAL MONTE CARLO MEANS ================\n")
print(
  main_table_means,
  row.names = FALSE
)

cat("\n================ FINAL PAPER TABLE: ALL METRICS ================\n")
print(
  paper_table_mean_sd,
  row.names = FALSE
)

cat(
  "\neta was saved as a descriptive diagnostic; ",
  "no single mixed-signal bias target is imposed.\n",
  sep = ""
)

cat("\nAll outputs were written to:\n")
cat(
  normalizePath(
    OUT_DIR,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)

cat("\n================ DONE ================\n")
