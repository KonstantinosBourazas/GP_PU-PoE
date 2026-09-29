## ============================================================================
## SIMULATION STUDY: SELF-CALIBRATING FULL PoE MODEL UNDER SAR
## ----------------------------------------------------------------------------
## TWO-EXPERT MODEL
## ----------------
## The fitted model is unchanged:
##
##   f = beta0 + X gamma_cov + Z_net gamma_net + g_cov + g_net,
##
##   g_cov ~ GP(0, K_RBF),
##   g_net ~ GP(0, K_heat),
##   Y_i | f_i ~ Bernoulli(plogis(f_i)).
##
## This is the ordinary Bernoulli-logistic two-expert PoE model. There is no
## PU parameter and no eta update.
##
## SAR LABEL MECHANISM
## -------------------
## The latent mixed-signal truth and every model, prior, MCMC setting, metric
## and output are retained from the SCAR script. The only substantive change is
## the observed-label mechanism:
##
##   - the 45 true positives remain fixed;
##   - exactly 15 are hidden without replacement in every replication;
##   - hiding weights are proportional to exp(-lambda * z_i), lambda = 2;
##   - z_i is the clipped standardized oracle covariate fraud-likeness score;
##   - the oracle score is used only to generate labels and is never supplied
##     to the fitted PoE model.
##
## The same SAR construction is used for the five global pilot data sets and
## the 100 production replications. Overall evaluation always uses the fixed
## latent truth T, so both observed and hidden positives have true label one.
##
## PRIOR CALIBRATION, MCMC AND EVALUATION
## --------------------------------------
## Geometry-specific prior calibration, the five global pilots, the exact dense
## two-expert GP sampler, the 15,000 production burn-in iterations, the 40,000
## posterior-sampling iterations thinned every 20, and all overall/hidden point
## and predictive-interval metrics are unchanged.
## ============================================================================

options(stringsAsFactors = FALSE)

## ============================================================================
## 0. PACKAGE CHECKS
## ============================================================================

required_namespaces <- c(
  "Matrix",
  "igraph",
  "RSpectra"
)

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
  library(igraph)
})

## ============================================================================
## 1. USER SETTINGS
## ============================================================================

SCRIPT_VERSION <- "full_poe_self_calibrating_priors_SAR_lambda2_v1_0"
EVALUATION_VERSION <- "node_predictive_labels_80_95_v2_0"

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
## Geometry-specific GP and structured-mean prior calibration
## --------------------------------------------------------------------------

## The active prior is recalculated separately for every pilot and production
## geometry. Unrounded values are always saved. To match the paper/application
## reporting convention, the active sampler parameters are rounded to two
## decimals after the exact geometry-specific calculation.
ROUND_ACTIVE_PRIORS <- TRUE
PRIOR_ROUND_DIGITS <- 2L

CORR_PCT_RANGE_ALL <- c(0.10, 0.50)
ELL_INTERVAL_PROBS <- c(0.025, 0.975)
SIGMA_INTERVAL_PROBS <- c(0.025, 0.975)
N_ELL_GRID_COARSE <- 60L
ROOT_TOL <- 1e-8

## Two experts split total latent residual variance equally.
TOTAL_SD_CENTER <- 2
COV_VARIANCE_SHARE <- 0.50
VARIANCE_SHARES <- c(
  cov = COV_VARIANCE_SHARE,
  network = 1 - COV_VARIANCE_SHARE
)
EFFECTIVE_SD_CENTER <- TOTAL_SD_CENTER * sqrt(VARIANCE_SHARES)
SIGMA_RANGE_FACTOR <- 2
EFFECTIVE_SD_LOW <- EFFECTIVE_SD_CENTER / SIGMA_RANGE_FACTOR
EFFECTIVE_SD_HIGH <- EFFECTIVE_SD_CENTER * SIGMA_RANGE_FACTOR

## Structured-mean half-t calibration, identical in form to the application.
HALFT_DF <- 3
MEAN_RMS_95_MULT <- 4.16
NU_TAU <- HALFT_DF
SAMPLE_TAU_GROUPS <- TRUE

## Fixed beta0 prior rule.
BETA0_PREV_RANGE_95 <- c(0.01, 0.10)

## Proposal and adaptation controls.
KAPPA_COV_INIT <- 1
KAPPA_NETWORK_INIT <- 1
PROPOSAL_COR_INIT <- list(
  cov = diag(2),
  network = diag(2)
)

ADAPT_BLOCK_CORRELATIONS_PILOTS <- TRUE
ADAPT_BLOCK_CORRELATIONS_PRODUCTION <- FALSE
COR_ADAPT_START <- 2000L
COR_ADAPT_INTERVAL <- 250L
COR_ADAPT_WINDOW <- 5000L
COR_SHRINKAGE <- 0.05
COR_MAX_ABS <- 0.95

EMP_COV_TAIL <- 5000L
JOINT_COR_SHRINKAGE <- 0.10
JOINT_COR_MAX_ABS <- 0.95
JOINT_THETA_SCALE_INIT_PILOT <- 1 / sqrt(2)
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
COV_JITTER <- 1e-6

JITTER <- 1e-8

## Common hard bounds used by the exact samplers and the calibration search.
MAX_LOG_SIGMA_COV <- log(5)
MAX_MEAN_SD_NETWORK <- 5
MAX_LOG_SIGMA_NETWORK <- log(1e3)
MIN_LOG_SIGMA <- log(1e-3)
LOGL_COV_BOUNDS <- c(log(1e-3), log(1e3))
LOGL_NETWORK_BOUNDS <- c(log(1e-4), log(1e2))
ELL_BOUNDS_COV <- exp(LOGL_COV_BOUNDS)
ELL_BOUNDS_NET <- exp(LOGL_NETWORK_BOUNDS)

## No spike component in the main simulation model.
PI_SPIKE <- 0
SPIKE_LOG_SIG_MEAN <- -6
SPIKE_LOG_SIG_SD <- 0.5

PILOT_PROGRESS_EVERY <- 5000L
PRODUCTION_BURNIN_PROGRESS_EVERY <- 1000L
PRODUCTION_SAMPLE_PROGRESS_EVERY <- 5000L

MOD_EIG_TOL <- 1e-8
MOD_MAX_Q <- Inf
MOD_K_MAX <- 20L

HIDDEN_TOP_K <- 15L
OVERALL_RECALL_K <- 45L
PREDICTIVE_INTERVAL_LEVELS <- c(0.80, 0.95)
PROBABILITY_EPS <- 1e-8

OUT_DIR <- file.path(
  getwd(),
  "SIMULATION_100_FULL_POE_SELF_CALIBRATING_PRIORS_SAR_LAMBDA_2"
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
  "PoE_global_pilot_calibration_5_datasets.rds"
)
GLOBAL_CALIBRATION_RDATA <- file.path(
  OUT_DIR,
  "PoE_global_pilot_calibration_5_datasets.RData"
)

RESUME_PILOTS <- TRUE
RESUME_PRODUCTION <- TRUE
FORCE_FRESH_PILOTS <- FALSE
FORCE_FRESH_PRODUCTION <- FALSE
FAIL_IF_PRODUCTION_FAILS <- TRUE
SAVE_PILOT_GEOMETRIES <- TRUE

## Reuse the SAR data geometries already seen by the competing methods
## when available. If absent, the same deterministic covariate, network,
## modularity and SAR seeds regenerate them internally.
USE_SAVED_COMPETING_GEOMETRIES <- TRUE
COMPETING_OUT_DIR <- file.path(
  getwd(),
  "SIMULATION_100_COMPETING_METHODS_MIXED_SIGNAL_SAR_LAMBDA_2"
)
COMPETING_GEOMETRY_DIR <- file.path(
  COMPETING_OUT_DIR,
  "generated_geometries"
)

## Deterministic independent MCMC streams.
PILOT_MCMC_SEED_BASE <- 710000000L
PRODUCTION_MCMC_SEED_BASE <- 810000000L
PRODUCTION_PREDICTIVE_SEED_BASE <- 910000000L

## ============================================================================
## 1A. FIXED LATENT TRUTH, MIXED-SIGNAL GEOMETRY AND SAR SETTINGS
## ============================================================================

N_NODES <- 400L
TRUE_ZERO_IDX <- 1:355
TRUE_POS_IDX <- 356:400

## These are the three positive signal groups used in the SCAR design.
## Their hidden/observed status is regenerated in every pilot and replication.
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

METHOD_NAME <- "PoE"
MODEL_KEY <- "full_poe"

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
## 1C. NETWORK DGP SETTINGS, NO SBM LAYER
## ============================================================================

ALPHA_BACKGROUND <- 2.0
ALPHA_FRAUD <- 2.0
ALPHA_LEGAL <- 2.0

D_BACKGROUND <- 2
D_FRAUD <- 2
D_LEGAL <- 2

TRIAD_CLOSURE_PROB <- 0.05
TRIAD_MAX_NEW_EDGES <- Inf

## ============================================================================
## 1D. REPRODUCIBLE SEED TABLES
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
  network_seed = as.integer(
    MASTER_SEED + 200000L + 1000L * seq_len(N_REP) + 29L
  ),
  modularity_seed = as.integer(
    MASTER_SEED + 300000L + 1000L * seq_len(N_REP) + 43L
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
  network_seed = as.integer(
    MASTER_SEED + 12000000L + 1000L * seq_len(N_PILOT) + 29L
  ),
  modularity_seed = as.integer(
    MASTER_SEED + 13000000L + 1000L * seq_len(N_PILOT) + 43L
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
  file.path(OUT_DIR, "PoE_production_seed_table.csv"),
  row.names = FALSE
)
write.csv(
  pilot_seed_table,
  file.path(OUT_DIR, "PoE_pilot_seed_table.csv"),
  row.names = FALSE
)

## ============================================================================
## 2. FIXED beta0 PRIOR AND FIXED COVARIATE TAU SCALE
## ============================================================================

round_prior <- function(x) {
  round(as.numeric(x), PRIOR_ROUND_DIGITS)
}

logit <- function(p) {
  log(p / (1 - p))
}

beta0_low <- logit(BETA0_PREV_RANGE_95[1L])
beta0_high <- logit(BETA0_PREV_RANGE_95[2L])

beta0_prior_for_sampler <- list(
  mean = round_prior(0.5 * (beta0_low + beta0_high)),
  sd = round_prior(
    (beta0_high - beta0_low) /
      (2 * qnorm(0.975))
  )
)

if (!is.finite(beta0_prior_for_sampler$mean) ||
    !is.finite(beta0_prior_for_sampler$sd) ||
    beta0_prior_for_sampler$sd <= 0) {
  stop("The internally constructed beta0 prior is invalid.")
}

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
  MEAN_RMS_95_MULT * unname(EFFECTIVE_SD_CENTER["cov"])
) / (
  HALF_T_Q95_UNIT * B_COV_FIXED
)

TAU_SCALE_COV_ACTIVE <- if (isTRUE(ROUND_ACTIVE_PRIORS)) {
  round_prior(TAU_SCALE_COV_UNROUNDED)
} else {
  TAU_SCALE_COV_UNROUNDED
}

if (!is.finite(TAU_SCALE_COV_ACTIVE) || TAU_SCALE_COV_ACTIVE <= 0) {
  stop("The fixed covariate half-t scale is invalid.")
}

cat("\n================ INTERNALLY CONSTRUCTED FIXED PRIOR ================\n")
cat(
  "beta0 ~ N(",
  beta0_prior_for_sampler$mean,
  ", ",
  beta0_prior_for_sampler$sd,
  "^2)\n",
  sep = ""
)
cat("Effective SD centers:\n")
print(EFFECTIVE_SD_CENTER)
cat("Fixed covariate mean dimension:", COVARIATE_MEAN_DIMENSION, "\n")
cat("Fixed covariate RMS factor     :", B_COV_FIXED, "\n")
cat("Fixed tau_cov half-t scale     :", TAU_SCALE_COV_ACTIVE, "\n")
cat("Active prior rounding          :", ROUND_ACTIVE_PRIORS, "\n")
cat("Active prior digits            :", PRIOR_ROUND_DIGITS, "\n")

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
  !identical(TRUE_POS_IDX, 356:400) ||
  !identical(which(TRUE_T == 1L), TRUE_POS_IDX)
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

if (!identical(
  sort(COVARIATE_SIGNAL_IDX),
  sort(c(
    356:360,
    366:380,
    391:400
  ))
)) {
  stop("COVARIATE_SIGNAL_IDX no longer matches the original SCAR geometry.")
}

if (!identical(
  sort(NETWORK_SIGNAL_IDX),
  sort(c(
    361:370,
    381:400
  ))
)) {
  stop("NETWORK_SIGNAL_IDX no longer matches the original SCAR geometry.")
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

if (PILOT_BURNIN <= 0L || PILOT_TAIL_LENGTH <= 0L ||
    PILOT_TAIL_LENGTH > PILOT_BURNIN) {
  stop("Invalid pilot warm-up or tail length.")
}
if (PRODUCTION_BURNIN <= 0L || PRODUCTION_SAMPLING <= 0L ||
    PRODUCTION_THIN <= 0L ||
    PRODUCTION_SAMPLING %% PRODUCTION_THIN != 0L) {
  stop("Invalid production MCMC lengths.")
}
if (!isTRUE(all.equal(PILOT_BLOCKED_FRAC, 1 / 3)) ||
    !isTRUE(all.equal(PRODUCTION_BLOCKED_FRAC, 1 / 3))) {
  stop("Pilot and production blocked fractions must both equal 1/3.")
}
if (length(PREDICTIVE_INTERVAL_LEVELS) != 2L ||
    !isTRUE(all.equal(PREDICTIVE_INTERVAL_LEVELS, c(0.80, 0.95))) ||
    any(PREDICTIVE_INTERVAL_LEVELS <= 0) ||
    any(PREDICTIVE_INTERVAL_LEVELS >= 1)) {
  stop("PREDICTIVE_INTERVAL_LEVELS must equal c(0.80, 0.95).")
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
    ELL_INTERVAL_PROBS[1L] <= 0 || ELL_INTERVAL_PROBS[2L] >= 1 ||
    SIGMA_INTERVAL_PROBS[1L] <= 0 || SIGMA_INTERVAL_PROBS[2L] >= 1) {
  stop("The ell and sigma interval probabilities are invalid.")
}
if (!isTRUE(all.equal(sum(VARIANCE_SHARES), 1)) ||
    any(VARIANCE_SHARES <= 0) ||
    !isTRUE(all.equal(
      sqrt(sum(EFFECTIVE_SD_CENTER^2)),
      TOTAL_SD_CENTER
    ))) {
  stop("The two-expert variance allocation is inconsistent.")
}
if (!identical(as.integer(NU_TAU), as.integer(HALFT_DF))) {
  stop("NU_TAU and HALFT_DF must be identical.")
}

if (isTRUE(FORCE_FRESH_PILOTS)) {
  unlink(PILOT_CHECKPOINT_DIR, recursive = TRUE, force = TRUE)
  dir.create(PILOT_CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(GLOBAL_CALIBRATION_FILE)) file.remove(GLOBAL_CALIBRATION_FILE)
  if (file.exists(GLOBAL_CALIBRATION_RDATA)) file.remove(GLOBAL_CALIBRATION_RDATA)
}
if (isTRUE(FORCE_FRESH_PRODUCTION)) {
  unlink(PRODUCTION_CHECKPOINT_DIR, recursive = TRUE, force = TRUE)
  dir.create(PRODUCTION_CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)
}

cat("\n================ FULL POE SAR SIMULATION DESIGN ================\n")
cat("Experts                  : covariate + one generated network\n")
cat("Pilot datasets           :", N_PILOT, "\n")
cat("Pilot warm-up            :", PILOT_BURNIN, "\n")
cat("Pilot pooled tail        :", PILOT_TAIL_LENGTH, "per pilot\n")
cat("Production replications  :", N_REP, "\n")
cat("Production burn-in       :", PRODUCTION_BURNIN, "\n")
cat("Production sampling      :", PRODUCTION_SAMPLING, "\n")
cat("Production thinning      :", PRODUCTION_THIN, "\n")
cat("Retained draws/run       :",
    PRODUCTION_SAMPLING / PRODUCTION_THIN, "\n")
cat("Exact GP matrix size     :", N_NODES, "x", N_NODES, "\n")
cat("Predictive interval levels:",
    paste(PREDICTIVE_INTERVAL_LEVELS, collapse = ", "), "\n")
cat("Prior calibration       : exact, geometry-specific in every run\n")
cat("GP effective SD centers :", paste(EFFECTIVE_SD_CENTER, collapse = ", "), "\n")
cat("tau half-t df           :", HALFT_DF, "\n")
cat("SAR mechanism           : covariate-driven, lambda =", LAMBDA_SAR, "\n")
cat("SAR z-score clipping    :", SAR_Z_CLIP, "\n")

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

logistic <- function(x) {
  plogis(x)
}

safe_logit <- function(p, eps = 1e-10) {
  qlogis(pmin(pmax(as.numeric(p), eps), 1 - eps))
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

sanitize_network_adjacency <- function(A) {
  A <- as.matrix(A)
  A[!is.finite(A)] <- 0
  A <- 1L * ((A + t(A)) > 0)
  diag(A) <- 0L
  A
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

sanitize_correlation_matrix <- function(
    R,
    max_abs = 0.95,
    shrinkage = 0,
    eig_floor = 1e-6
) {
  R <- as.matrix(R)
  if (nrow(R) != ncol(R) || nrow(R) < 1L) {
    stop("R must be a non-empty square matrix.")
  }
  p <- nrow(R)
  R[!is.finite(R)] <- 0
  R <- symmetrize(R)
  diag(R) <- 1
  if (p > 1L) {
    off <- row(R) != col(R)
    R[off] <- pmax(-abs(max_abs), pmin(abs(max_abs), R[off]))
  }
  shrinkage <- min(max(as.numeric(shrinkage), 0), 1)
  R <- (1 - shrinkage) * R + shrinkage * diag(p)
  ee <- eigen(symmetrize(R), symmetric = TRUE)
  vals <- pmax(Re(ee$values), eig_floor)
  Rpd <- Re(ee$vectors %*% (vals * t(ee$vectors)))
  d <- sqrt(pmax(diag(Rpd), eig_floor))
  Rpd <- Rpd / outer(d, d)
  Rpd <- symmetrize(Rpd)
  diag(Rpd) <- 1
  Rpd
}

estimate_correlation_matrix <- function(
    M,
    fallback = NULL,
    shrinkage = 0.10,
    max_abs = 0.95,
    min_n = 100L
) {
  M <- as.matrix(M)
  M <- M[complete.cases(M), , drop = FALSE]
  p <- ncol(M)
  if (is.null(fallback)) fallback <- diag(p)
  if (nrow(M) < as.integer(min_n) || p < 1L) {
    return(sanitize_correlation_matrix(
      fallback,
      max_abs = max_abs,
      shrinkage = shrinkage
    ))
  }
  sdv <- apply(M, 2L, sd)
  if (any(!is.finite(sdv)) || any(sdv <= 0)) {
    return(sanitize_correlation_matrix(
      fallback,
      max_abs = max_abs,
      shrinkage = shrinkage
    ))
  }
  R <- suppressWarnings(cor(M))
  sanitize_correlation_matrix(
    R,
    max_abs = max_abs,
    shrinkage = shrinkage
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

simulate_binary_label_draws <- function(probability_matrix, seed) {
  probability_matrix <- as.matrix(probability_matrix)
  if (nrow(probability_matrix) < 1L || ncol(probability_matrix) < 1L ||
      any(!is.finite(probability_matrix)) ||
      any(probability_matrix < 0) ||
      any(probability_matrix > 1)) {
    stop("The predictive probability matrix is invalid.")
  }
  
  with_local_seed(seed, {
    draws <- matrix(
      as.integer(
        runif(length(probability_matrix)) <
          as.vector(probability_matrix)
      ),
      nrow = nrow(probability_matrix),
      ncol = ncol(probability_matrix)
    )
    dimnames(draws) <- dimnames(probability_matrix)
    draws
  })
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

overall_brier <- function(probabilities, truth) {
  mean((as.numeric(probabilities) - as.numeric(truth))^2)
}

interval_score_vector <- function(lower, upper, truth, alpha) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  truth <- as.numeric(truth)
  alpha <- as.numeric(alpha)
  
  if (length(lower) != length(upper) ||
      length(lower) != length(truth) ||
      length(alpha) != 1L ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1) {
    stop("Invalid inputs to interval_score_vector.")
  }
  
  (upper - lower) +
    (2 / alpha) * (lower - truth) * (truth < lower) +
    (2 / alpha) * (truth - upper) * (truth > upper)
}

predictive_set_metrics <- function(
    predictive_T_draws,
    truth,
    evaluation_idx,
    level
) {
  predictive_T_draws <- as.matrix(predictive_T_draws)
  truth <- as.integer(truth)
  evaluation_idx <- as.integer(evaluation_idx)
  
  if (ncol(predictive_T_draws) != length(truth)) {
    stop("predictive_T_draws and truth have incompatible dimensions.")
  }
  if (any(!predictive_T_draws %in% c(0L, 1L))) {
    stop("predictive_T_draws must be binary.")
  }
  if (length(evaluation_idx) < 1L ||
      any(evaluation_idx < 1L) ||
      any(evaluation_idx > length(truth))) {
    stop("evaluation_idx is invalid.")
  }
  if (!is.finite(level) || level <= 0 || level >= 1) {
    stop("level must lie strictly between zero and one.")
  }
  
  alpha <- 1 - level
  draws_subset <- predictive_T_draws[
    ,
    evaluation_idx,
    drop = FALSE
  ]
  truth_subset <- truth[evaluation_idx]
  
  lower <- apply(
    draws_subset,
    2L,
    stats::quantile,
    probs = alpha / 2,
    type = 1,
    names = FALSE,
    na.rm = TRUE
  )
  upper <- apply(
    draws_subset,
    2L,
    stats::quantile,
    probs = 1 - alpha / 2,
    type = 1,
    names = FALSE,
    na.rm = TRUE
  )
  
  covered <- as.numeric(
    lower <= truth_subset &
      truth_subset <= upper
  )
  scores <- interval_score_vector(
    lower = lower,
    upper = upper,
    truth = truth_subset,
    alpha = alpha
  )
  
  list(
    level = level,
    alpha = alpha,
    evaluation_idx = evaluation_idx,
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    covered = covered,
    interval_score = as.numeric(scores),
    coverage = mean(covered),
    mean_interval_score = mean(scores),
    mean_width = mean(upper - lower)
  )
}

summary_stats <- function(x) {
  x <- as.numeric(x)
  n_ok <- sum(is.finite(x))
  c(
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    mc_se = stats::sd(x, na.rm = TRUE) / sqrt(n_ok),
    q025 = unname(stats::quantile(x, 0.025, na.rm = TRUE)),
    median = stats::median(x, na.rm = TRUE),
    q975 = unname(stats::quantile(x, 0.975, na.rm = TRUE))
  )
}

file_md5_or_na <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path)[1L])
}

normalize_named_positive <- function(x, group_names, arg_name) {
  x_names <- names(x)
  x_num <- as.numeric(x)
  if (length(x_num) != length(group_names)) {
    stop(arg_name, " has the wrong length.")
  }
  if (!is.null(x_names)) {
    if (!all(group_names %in% x_names)) {
      stop(arg_name, " misses required group names.")
    }
    out <- as.numeric(x[group_names])
  } else {
    out <- x_num
  }
  names(out) <- group_names
  if (any(!is.finite(out)) || any(out <= 0)) {
    stop(arg_name, " must contain positive finite values.")
  }
  out
}

median_covariance_matrix <- function(covariance_list, dimension) {
  if (length(covariance_list) < 1L) {
    stop("At least one covariance matrix is required.")
  }
  mats <- lapply(covariance_list, function(S) {
    S <- as.matrix(S)
    if (!all(dim(S) == c(dimension, dimension)) || any(!is.finite(S))) {
      stop("All pilot covariance matrices have to be finite and conformable.")
    }
    S
  })
  out <- matrix(0, dimension, dimension)
  for (ii in seq_len(dimension)) {
    for (jj in ii:dimension) {
      value <- median(vapply(mats, function(S) S[ii, jj], numeric(1)))
      out[ii, jj] <- value
      out[jj, ii] <- value
    }
  }
  sanitize_covariance_matrix(out + COV_JITTER * diag(dimension))
}
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
## 4. COVARIATE DATA-GENERATING PROCESS
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
## 5. NETWORK DATA-GENERATING PROCESS
## ============================================================================

edge_count <- function(A) {
  sum(A) / 2
}

add_edges <- function(
    A,
    pairs
) {
  if (
    is.null(pairs) ||
    length(pairs) == 0L ||
    nrow(as.matrix(pairs)) == 0L
  ) {
    return(A)
  }
  
  pairs <- as.matrix(pairs)
  
  A[
    cbind(
      pairs[, 1L],
      pairs[, 2L]
    )
  ] <- 1L
  
  A[
    cbind(
      pairs[, 2L],
      pairs[, 1L]
    )
  ] <- 1L
  
  diag(A) <- 0L
  
  A
}

crp_partition <- function(
    nodes,
    alpha
) {
  nodes <- as.integer(nodes)
  
  if (length(nodes) == 0L) {
    return(
      list()
    )
  }
  
  if (
    !is.finite(alpha) ||
    alpha <= 0
  ) {
    stop(
      "alpha must be positive."
    )
  }
  
  order_in <- sample(
    nodes,
    length(nodes),
    replace = FALSE
  )
  
  groups <- list(
    order_in[1L]
  )
  
  if (length(order_in) >= 2L) {
    for (
      ii in 2:length(order_in)
    ) {
      group_sizes <- lengths(
        groups
      )
      
      choice <- sample.int(
        length(group_sizes) + 1L,
        size = 1L,
        prob = c(
          group_sizes,
          alpha
        )
      )
      
      if (
        choice ==
        length(group_sizes) + 1L
      ) {
        groups[[length(groups) + 1L]] <- order_in[ii]
      } else {
        groups[[choice]] <- c(
          groups[[choice]],
          order_in[ii]
        )
      }
    }
  }
  
  names(groups) <- paste0(
    "G",
    seq_along(groups)
  )
  
  groups
}

build_crp_edge_layer <- function(
    n_nodes,
    groups,
    target_degree
) {
  A <- matrix(
    0L,
    n_nodes,
    n_nodes
  )
  
  for (
    group_index in seq_along(groups)
  ) {
    members <- as.integer(
      groups[[group_index]]
    )
    
    group_size <- length(
      members
    )
    
    if (group_size < 2L) {
      next
    }
    
    edge_probability <- min(
      1,
      target_degree /
        (
          group_size - 1
        )
    )
    
    candidates <- t(
      utils::combn(
        members,
        2L
      )
    )
    
    keep <- runif(
      nrow(candidates)
    ) < edge_probability
    
    if (any(keep)) {
      A <- add_edges(
        A,
        candidates[
          keep,
          ,
          drop = FALSE
        ]
      )
    }
  }
  
  sanitize_network_adjacency(
    A
  )
}

apply_triadic_closure <- function(
    A_current,
    closure_prob = TRIAD_CLOSURE_PROB,
    max_new_edges = TRIAD_MAX_NEW_EDGES
) {
  A_current <- sanitize_network_adjacency(
    A_current
  )
  
  n <- nrow(
    A_current
  )
  
  common_neighbors <- A_current %*%
    A_current
  
  candidate_mask <- upper.tri(
    A_current,
    diag = FALSE
  ) &
    A_current == 0L &
    common_neighbors > 0
  
  candidate_pairs <- which(
    candidate_mask,
    arr.ind = TRUE
  )
  
  if (
    nrow(candidate_pairs) == 0L ||
    closure_prob <= 0
  ) {
    return(
      list(
        A = matrix(
          0L,
          n,
          n
        ),
        n_candidates = 0L,
        n_edges = 0L
      )
    )
  }
  
  common_neighbor_count <- common_neighbors[
    cbind(
      candidate_pairs[, 1L],
      candidate_pairs[, 2L]
    )
  ]
  
  pair_probability <- 1 -
    (
      1 - closure_prob
    )^common_neighbor_count
  
  keep <- runif(
    length(pair_probability)
  ) < pair_probability
  
  chosen_pairs <- candidate_pairs[
    keep,
    ,
    drop = FALSE
  ]
  
  if (
    is.finite(max_new_edges) &&
    nrow(chosen_pairs) >
    max_new_edges
  ) {
    chosen_pairs <- chosen_pairs[
      sample.int(
        nrow(chosen_pairs),
        max_new_edges
      ),
      ,
      drop = FALSE
    ]
  }
  
  A_closure <- matrix(
    0L,
    n,
    n
  )
  
  A_closure <- add_edges(
    A_closure,
    chosen_pairs
  )
  
  list(
    A = sanitize_network_adjacency(
      A_closure
    ),
    n_candidates = nrow(
      candidate_pairs
    ),
    n_edges = edge_count(
      A_closure
    )
  )
}

generate_network_set <- function(
    seed,
    signal_idx = NETWORK_SIGNAL_IDX
) {
  set.seed(seed)
  
  signal_idx <- sort(
    unique(
      as.integer(signal_idx)
    )
  )
  
  if (length(signal_idx) != 30L) {
    stop(
      "The network DGP requires exactly 30 signal nodes."
    )
  }
  
  background_idx <- setdiff(
    seq_len(N_NODES),
    signal_idx
  )
  
  background_groups <- crp_partition(
    nodes = seq_len(
      N_NODES
    ),
    alpha = ALPHA_BACKGROUND
  )
  
  A_background <- build_crp_edge_layer(
    n_nodes = N_NODES,
    groups = background_groups,
    target_degree = D_BACKGROUND
  )
  
  signal_groups <- crp_partition(
    nodes = signal_idx,
    alpha = ALPHA_FRAUD
  )
  
  A_signal_raw <- build_crp_edge_layer(
    n_nodes = N_NODES,
    groups = signal_groups,
    target_degree = D_FRAUD
  )
  
  A_after_signal <- sanitize_network_adjacency(
    A_background +
      A_signal_raw
  )
  
  legal_groups <- crp_partition(
    nodes = background_idx,
    alpha = ALPHA_LEGAL
  )
  
  A_legal_raw <- build_crp_edge_layer(
    n_nodes = N_NODES,
    groups = legal_groups,
    target_degree = D_LEGAL
  )
  
  A_local <- sanitize_network_adjacency(
    A_after_signal +
      A_legal_raw
  )
  
  local_edges <- edge_count(
    A_local
  )
  
  ## The final network before closure is the union of the three CRP layers.
  A_preclosure <- A_local
  
  closure <- apply_triadic_closure(
    A_current = A_preclosure,
    closure_prob = TRIAD_CLOSURE_PROB,
    max_new_edges = TRIAD_MAX_NEW_EDGES
  )
  
  A_final <- sanitize_network_adjacency(
    A_preclosure +
      closure$A
  )
  
  degree <- rowSums(
    A_final
  )
  
  graph_final <- igraph::graph_from_adjacency_matrix(
    A_final,
    mode = "undirected",
    diag = FALSE
  )
  
  component_info <- igraph::components(
    graph_final
  )
  
  list(
    seed = seed,
    A = A_final,
    signal_idx = signal_idx,
    background_idx = background_idx,
    background_groups = background_groups,
    signal_groups = signal_groups,
    legal_groups = legal_groups,
    degree = degree,
    summary = data.frame(
      n_nodes = N_NODES,
      n_network_signal = length(
        signal_idx
      ),
      n_network_background = length(
        background_idx
      ),
      n_background_groups = length(
        background_groups
      ),
      n_signal_groups = length(
        signal_groups
      ),
      n_legal_groups = length(
        legal_groups
      ),
      edges_local = local_edges,
      edges_closure = edge_count(
        closure$A
      ),
      edges_final = edge_count(
        A_final
      ),
      mean_degree = mean(
        degree
      ),
      mean_degree_network_background = mean(
        degree[
          background_idx
        ]
      ),
      mean_degree_network_signal = mean(
        degree[
          signal_idx
        ]
      ),
      n_components = component_info$no,
      largest_component = max(
        component_info$csize
      ),
      stringsAsFactors = FALSE
    )
  )
}

## ============================================================================
## 6. LABEL-BLIND MODULARITY FEATURES
## ============================================================================

choose_q_modularity_eigengap <- function(
    values,
    tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q
) {
  positive_values <- sort(
    as.numeric(
      values[
        values > tol
      ]
    ),
    decreasing = TRUE
  )
  
  n_positive <- length(
    positive_values
  )
  
  if (n_positive == 0L) {
    return(
      list(
        q = 1L,
        positive_values = numeric(0),
        gaps = numeric(0),
        max_gap = NA_real_
      )
    )
  }
  
  if (n_positive == 1L) {
    return(
      list(
        q = 1L,
        positive_values = positive_values,
        gaps = numeric(0),
        max_gap = NA_real_
      )
    )
  }
  
  gaps <- positive_values[
    -n_positive
  ] -
    positive_values[
      -1L
    ]
  
  q <- which.max(
    gaps
  )
  
  if (is.finite(max_q)) {
    q <- min(
      q,
      as.integer(max_q)
    )
  }
  
  q <- max(
    1L,
    min(
      q,
      n_positive
    )
  )
  
  list(
    q = as.integer(q),
    positive_values = positive_values,
    gaps = gaps,
    max_gap = gaps[
      which.max(gaps)
    ]
  )
}

orient_eigenvectors_deterministically <- function(V) {
  V <- as.matrix(V)
  
  for (
    column_index in seq_len(
      ncol(V)
    )
  ) {
    largest_loading_index <- which.max(
      abs(
        V[
          ,
          column_index
        ]
      )
    )
    
    if (
      V[
        largest_loading_index,
        column_index
      ] < 0
    ) {
      V[
        ,
        column_index
      ] <- -V[
        ,
        column_index
      ]
    }
  }
  
  V
}

modularity_features_single_network <- function(
    A,
    seed,
    eig_tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q,
    k_max = MOD_K_MAX
) {
  set.seed(seed)
  
  A <- sanitize_network_adjacency(
    A
  )
  
  n <- nrow(
    A
  )
  
  degree <- rowSums(
    A
  )
  
  n_edges <- sum(
    degree
  ) / 2
  
  if (
    !is.finite(n_edges) ||
    n_edges <= 0
  ) {
    stop(
      "Cannot construct modularity features from an edgeless network."
    )
  }
  
  modularity_operator <- function(
    x,
    args
  ) {
    as.vector(
      A %*% x -
        degree *
        sum(
          degree * x
        ) /
        (
          2 * n_edges
        )
    )
  }
  
  k_use <- min(
    as.integer(k_max),
    n - 2L
  )
  
  eig <- RSpectra::eigs_sym(
    A = modularity_operator,
    k = k_use,
    n = n,
    which = "LA",
    opts = list(
      retvec = TRUE
    )
  )
  
  eigenvalues <- Re(
    eig$values
  )
  
  eigenvectors <- Re(
    eig$vectors
  )
  
  ordering <- order(
    eigenvalues,
    decreasing = TRUE
  )
  
  eigenvalues <- eigenvalues[
    ordering
  ]
  
  eigenvectors <- eigenvectors[
    ,
    ordering,
    drop = FALSE
  ]
  
  eigenvectors <- orient_eigenvectors_deterministically(
    eigenvectors
  )
  
  positive_index <- which(
    eigenvalues > eig_tol
  )
  
  if (length(positive_index) == 0L) {
    keep <- 1L
  } else {
    gap_information <- choose_q_modularity_eigengap(
      values = eigenvalues,
      tol = eig_tol,
      max_q = max_q
    )
    
    keep <- positive_index[
      seq_len(
        gap_information$q
      )
    ]
  }
  
  Z_raw <- eigenvectors[
    ,
    keep,
    drop = FALSE
  ]
  
  colnames(
    Z_raw
  ) <- paste0(
    "network_mod",
    seq_len(
      ncol(Z_raw)
    )
  )
  
  standardized <- standardize_columns(
    Z_raw
  )
  
  Z <- standardized$scaled
  
  colnames(Z) <- colnames(
    Z_raw
  )
  
  list(
    Z = Z,
    raw = Z_raw,
    q = ncol(Z),
    selected_eigenvalues = eigenvalues[
      keep
    ],
    all_eigenvalues = eigenvalues,
    center = standardized$center,
    scale = standardized$scale
  )
}

## ============================================================================
## 7. EXACT NETWORK-GP GEOMETRY
## ============================================================================

prepare_network_components_for_blocks <- function(
    A,
    network_name = "network"
) {
  A <- sanitize_network_adjacency(A)
  n <- nrow(A)
  
  g <- igraph::graph_from_adjacency_matrix(
    A,
    mode = "undirected",
    diag = FALSE
  )
  
  comp <- igraph::components(g)$membership
  comp_ids <- sort(unique(comp))
  components <- vector("list", length(comp_ids))
  
  for (cc_pos in seq_along(comp_ids)) {
    cc <- comp_ids[cc_pos]
    idx <- which(comp == cc)
    s <- length(idx)
    
    object <- list(idx = idx, size = s)
    
    if (s == 1L) {
      object$type <- "isolate"
      components[[cc_pos]] <- object
      next
    }
    
    A_sub <- A[idx, idx, drop = FALSE]
    edge_total <- sum(A_sub) / 2
    is_clique <- abs(edge_total - s * (s - 1) / 2) < 1e-8
    
    if (is_clique) {
      object$type <- "clique"
      object$lambda_orth <- s / (s - 1)
      components[[cc_pos]] <- object
      next
    }
    
    degree <- rowSums(A_sub)
    inverse_sqrt_degree <- rep(0, s)
    positive_degree <- which(degree > 0)
    inverse_sqrt_degree[positive_degree] <-
      1 / sqrt(degree[positive_degree])
    
    D_inverse <- diag(inverse_sqrt_degree, s)
    L_normalized <- diag(1, s) -
      D_inverse %*% A_sub %*% D_inverse
    L_normalized <- symmetrize(L_normalized)
    
    eig <- eigen(L_normalized, symmetric = TRUE)
    
    object$type <- "general"
    object$U <- Re(eig$vectors)
    object$U2 <- object$U^2
    object$lambda <- pmax(Re(eig$values), 0)
    components[[cc_pos]] <- object
  }
  
  component_sizes <- vapply(components, `[[`, integer(1), "size")
  component_types <- vapply(components, `[[`, character(1), "type")
  
  list(
    name = network_name,
    A = A,
    components = components,
    component_membership = comp,
    component_summary = data.frame(
      network = network_name,
      n_components = length(components),
      max_component_size = max(component_sizes),
      n_isolates = sum(component_types == "isolate"),
      n_cliques = sum(component_types == "clique"),
      n_general = sum(component_types == "general"),
      stringsAsFactors = FALSE
    )
  )
}

network_kernel_from_geom <- function(geom, ell, sigma) {
  n <- nrow(geom$A)
  K <- matrix(0, n, n)
  sigma2 <- sigma^2
  
  for (component in geom$components) {
    idx <- component$idx
    s <- component$size
    
    if (component$type == "isolate") {
      q <- exp(-1 / (2 * ell^2))
      K[idx, idx] <- K[idx, idx] + sigma2 * q
      next
    }
    
    if (component$type == "clique") {
      q <- exp(-component$lambda_orth / (2 * ell^2))
      off_diagonal <- sigma2 * (1 - q) / s
      K[idx, idx] <- K[idx, idx] + off_diagonal
      K[cbind(idx, idx)] <- K[cbind(idx, idx)] + sigma2 * q
      next
    }
    
    gamma <- sigma2 * exp(-component$lambda / (2 * ell^2))
    U_scaled <- sweep(component$U, 2L, sqrt(gamma), "*")
    K_component <- symmetrize(tcrossprod(U_scaled))
    K[idx, idx] <- K[idx, idx] + K_component
  }
  
  symmetrize(K)
}

## ============================================================================
## 8. EXACT GEOMETRY-SPECIFIC PRIOR CALIBRATION
## ============================================================================

fit_lognormal_interval <- function(
    q_low,
    q_high,
    probs = c(0.025, 0.975)
) {
  q_low <- as.numeric(q_low)
  q_high <- as.numeric(q_high)
  
  if (length(q_low) != 1L || length(q_high) != 1L ||
      !is.finite(q_low) || !is.finite(q_high) ||
      q_low <= 0 || q_high <= 0) {
    stop("q_low and q_high must be positive finite scalars.")
  }
  
  qlo <- min(q_low, q_high)
  qhi <- max(q_low, q_high)
  p_low <- probs[1L]
  p_high <- probs[2L]
  
  if (!is.finite(p_low) || !is.finite(p_high) ||
      p_low <= 0 || p_high >= 1 || p_low >= p_high) {
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
    if (length(ell) != 1L || !is.finite(ell) || ell <= 0) {
      return(list(mean_corr = NA_real_, mean_sd = NA_real_))
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
    if (length(ell) != 1L || !is.finite(ell) || ell <= 0) {
      return(list(mean_corr = NA_real_, mean_sd = NA_real_))
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
      
      if (is.finite(off_corr_sum) && off_corr_sum < 0 &&
          abs(off_corr_sum) < 1e-8) {
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
    is.finite(curve$ell) & is.finite(curve$mean_corr),
    ,
    drop = FALSE
  ]
  curve <- curve[order(curve$ell), , drop = FALSE]
  
  corr_range <- range(curve$mean_corr, na.rm = TRUE)
  target_used <- min(max(target_corr, corr_range[1L]), corr_range[2L])
  difference <- curve$mean_corr - target_used
  exact_idx <- which(abs(difference) <= 1e-12)
  
  if (length(exact_idx) > 0L) {
    ell_hat <- curve$ell[exact_idx[1L]]
  } else {
    crossing_idx <- which(
      difference[-nrow(curve)] * difference[-1L] <= 0
    )
    root_result <- NULL
    
    if (length(crossing_idx) > 0L) {
      crossing_score <- abs(difference[crossing_idx]) +
        abs(difference[crossing_idx + 1L])
      kk <- crossing_idx[which.min(crossing_score)]
      log_interval <- log(c(curve$ell[kk], curve$ell[kk + 1L]))
      
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
    stop("Failed to invert the correlation curve for ", expert_name, ".")
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
  
  ## Never allow normal rounding to collapse a positive prior SD to zero.
  if (output$log_ell_sd <= 0) {
    output$log_ell_sd <- hyperprior_unrounded$log_ell_sd
  }
  if (output$log_sig_sd <= 0) {
    output$log_sig_sd <- hyperprior_unrounded$log_sig_sd
  }
  
  if (any(!is.finite(unlist(output))) ||
      output$log_ell_sd <= 0 || output$log_sig_sd <= 0) {
    stop("The active geometry-specific hyperprior is invalid.")
  }
  
  output
}

calibrate_geometry_priors <- function(
    X_cov,
    Z_network,
    network_geometry,
    calibration_label
) {
  calibration_start <- Sys.time()
  X_cov <- as.matrix(X_cov)
  Z_network <- as.matrix(Z_network)
  
  if (nrow(X_cov) != N_NODES || ncol(X_cov) != 10L ||
      nrow(Z_network) != N_NODES || ncol(Z_network) < 1L) {
    stop("Geometry-specific prior calibration received invalid designs.")
  }
  
  covariate_metric <- make_covariate_metric_fun(X_cov)
  network_metric <- network_metric_fun_factory(network_geometry)
  
  covariate_calibration <- calibrate_one_expert_prior(
    expert_name = "cov",
    metric_fun = covariate_metric,
    ell_bounds = ELL_BOUNDS_COV,
    effective_sd_center = EFFECTIVE_SD_CENTER["cov"]
  )
  network_calibration <- calibrate_one_expert_prior(
    expert_name = "network",
    metric_fun = network_metric,
    ell_bounds = ELL_BOUNDS_NET,
    effective_sd_center = EFFECTIVE_SD_CENTER["network"]
  )
  
  hyperpriors_unrounded <- list(
    expert0_cov = covariate_calibration$hyperprior_unrounded,
    expert1_network = network_calibration$hyperprior_unrounded
  )
  active_hyperpriors <- lapply(
    hyperpriors_unrounded,
    activate_hyperprior
  )
  
  b_cov_actual <- sqrt(sum(X_cov^2) / nrow(X_cov))
  b_network <- sqrt(sum(Z_network^2) / nrow(Z_network))
  b_network_from_dimension <- sqrt(
    ncol(Z_network) *
      (nrow(Z_network) - 1) /
      nrow(Z_network)
  )
  
  if (!is.finite(b_cov_actual) || !is.finite(b_network) ||
      b_cov_actual <= 0 || b_network <= 0) {
    stop("Invalid structured-mean RMS factors.")
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
  if (!isTRUE(all.equal(
    b_network,
    b_network_from_dimension,
    tolerance = 1e-8
  ))) {
    stop(
      "The modularity design is not standardized as expected: b_network = ",
      signif(b_network, 10),
      ", expected from q = ",
      signif(b_network_from_dimension, 10),
      "."
    )
  }
  
  tau_scale_network_unrounded <- (
    MEAN_RMS_95_MULT * unname(EFFECTIVE_SD_CENTER["network"])
  ) / (
    HALF_T_Q95_UNIT * b_network
  )
  
  tau_scale_unrounded <- c(
    cov = TAU_SCALE_COV_UNROUNDED,
    network = tau_scale_network_unrounded
  )
  tau_scale_active <- if (isTRUE(ROUND_ACTIVE_PRIORS)) {
    round_prior(tau_scale_unrounded)
  } else {
    tau_scale_unrounded
  }
  names(tau_scale_active) <- names(tau_scale_unrounded)
  
  if (any(!is.finite(tau_scale_active)) ||
      any(tau_scale_active <= 0)) {
    stop("The geometry-specific half-t scales are invalid.")
  }
  
  tau2_default <- tau_scale_active^2
  
  theta_prior_mean <- c(
    log_ell_cov = active_hyperpriors$expert0_cov$log_ell_mean,
    log_sigma_cov = active_hyperpriors$expert0_cov$log_sig_mean,
    log_ell_network = active_hyperpriors$expert1_network$log_ell_mean,
    log_sigma_network = active_hyperpriors$expert1_network$log_sig_mean
  )
  theta_prior_sd <- c(
    log_ell_cov = active_hyperpriors$expert0_cov$log_ell_sd,
    log_sigma_cov = active_hyperpriors$expert0_cov$log_sig_sd,
    log_ell_network = active_hyperpriors$expert1_network$log_ell_sd,
    log_sigma_network = active_hyperpriors$expert1_network$log_sig_sd
  )
  
  summary_table <- rbind(
    covariate_calibration$row,
    network_calibration$row
  )
  summary_table$n_mean_features <- c(ncol(X_cov), ncol(Z_network))
  summary_table$Z_rms_factor <- c(b_cov_actual, b_network)
  summary_table$Z_rms_factor_from_dimension <- c(
    B_COV_FIXED,
    b_network_from_dimension
  )
  summary_table$half_t_df <- HALFT_DF
  summary_table$mean_rms_95_mult <- MEAN_RMS_95_MULT
  summary_table$tau_scale_unrounded <- unname(
    tau_scale_unrounded[summary_table$expert]
  )
  summary_table$tau_scale_active <- unname(
    tau_scale_active[summary_table$expert]
  )
  summary_table$tau2_default <- unname(
    tau2_default[summary_table$expert]
  )
  summary_table$log_ell_mean_active <- c(
    active_hyperpriors$expert0_cov$log_ell_mean,
    active_hyperpriors$expert1_network$log_ell_mean
  )
  summary_table$log_ell_sd_active <- c(
    active_hyperpriors$expert0_cov$log_ell_sd,
    active_hyperpriors$expert1_network$log_ell_sd
  )
  summary_table$log_sig_mean_active <- c(
    active_hyperpriors$expert0_cov$log_sig_mean,
    active_hyperpriors$expert1_network$log_sig_mean
  )
  summary_table$log_sig_sd_active <- c(
    active_hyperpriors$expert0_cov$log_sig_sd,
    active_hyperpriors$expert1_network$log_sig_sd
  )
  summary_table$calibration_label <- calibration_label
  summary_table$active_prior_rounding <- ROUND_ACTIVE_PRIORS
  summary_table$active_prior_digits <- PRIOR_ROUND_DIGITS
  
  runtime_sec <- as.numeric(difftime(
    Sys.time(),
    calibration_start,
    units = "secs"
  ))
  summary_table$prior_calibration_runtime_sec <- runtime_sec
  
  list(
    calibration_label = calibration_label,
    active_hyperpriors = active_hyperpriors,
    hyperpriors_unrounded = hyperpriors_unrounded,
    tau_scale = tau_scale_active,
    tau_scale_unrounded = tau_scale_unrounded,
    tau2_default = tau2_default,
    theta_prior_mean = theta_prior_mean,
    theta_prior_sd = theta_prior_sd,
    summary = summary_table,
    runtime_sec = runtime_sec
  )
}

## ============================================================================
## 9. BERNOULLI-LOGISTIC LIKELIHOOD
## ============================================================================

bernoulli_loglik_vec <- function(f, y) {
  f <- as.numeric(f)
  y <- as.numeric(y)
  softplus_f <- pmax(f, 0) + log1p(exp(-abs(f)))
  sum(y * f - softplus_f)
}

bernoulli_gradloglik_vec <- function(f, y) {
  as.numeric(y) - plogis(as.numeric(f))
}

## ============================================================================
## 10. EXACT TWO-EXPERT FULL PoE CHAIN
## ============================================================================

run_full_poe_exact_chain <- function(
    Y,
    X_cov,
    A_network,
    Z_network,
    network_geometry,
    hyperpriors,
    m_beta,
    s_beta,
    tau2_default,
    scale_tau,
    n_burnin,
    n_sample = 0L,
    thin = 1L,
    blocked_frac = 1 / 3,
    min_blocked_warmup = 5000L,
    init_state = NULL,
    delta_init = 3.5,
    kappa_cov_init = 1,
    kappa_network_init = 1,
    proposal_cor_init = NULL,
    adapt_block_cor = TRUE,
    cor_adapt_start = 2000L,
    cor_adapt_interval = 250L,
    cor_adapt_window = 5000L,
    cor_shrinkage = 0.05,
    cor_max_abs = 0.95,
    joint_cov_base_init = NULL,
    joint_theta_scale_init = 1 / sqrt(2),
    joint_cor_shrinkage = 0.10,
    joint_cor_max_abs = 0.95,
    f_refresh_after_joint = TRUE,
    nu_tau = 3,
    sample_tau_groups = TRUE,
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
    max_mean_sd_network = 5,
    max_log_sigma_network = log(1e3),
    min_log_sigma = log(1e-3),
    logl_cov_bounds = c(log(1e-3), log(1e3)),
    logl_network_bounds = c(log(1e-4), log(1e2)),
    pi_spike = 0,
    spike_log_sig_mean = -6,
    spike_log_sig_sd = 0.5,
    jitter = 1e-8,
    state_tail_length = 0L,
    predictive_label_seed = NULL,
    phase_label = "PoE",
    burnin_progress_every = 5000L,
    sample_progress_every = 5000L,
    verbose = TRUE
) {
  Y <- as.numeric(Y)
  X_cov <- as.matrix(X_cov)
  A_network <- sanitize_network_adjacency(A_network)
  Z_network <- as.matrix(Z_network)
  
  n <- length(Y)
  p_cov <- ncol(X_cov)
  q_network <- ncol(Z_network)
  
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
  if (nrow(X_cov) != n || nrow(Z_network) != n ||
      !all(dim(A_network) == c(n, n))) {
    stop("The data objects have incompatible dimensions.")
  }
  if (p_cov != 10L) {
    stop("The full PoE simulation expects exactly ten covariate columns.")
  }
  if (q_network < 1L) {
    stop("The network mean design must contain at least one modularity feature.")
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
  
  group_names <- c("cov", "network")
  theta_names <- c(
    "log_ell_cov",
    "log_sigma_cov",
    "log_ell_network",
    "log_sigma_network"
  )
  
  scale_tau <- normalize_named_positive(
    scale_tau,
    group_names,
    "scale_tau"
  )
  tau2_default <- normalize_named_positive(
    tau2_default,
    group_names,
    "tau2_default"
  )
  
  if (is.null(colnames(X_cov))) {
    colnames(X_cov) <- paste0("X", seq_len(p_cov))
  }
  if (is.null(colnames(Z_network))) {
    colnames(Z_network) <- paste0("network_mod", seq_len(q_network))
  }
  colnames(X_cov) <- make.unique(colnames(X_cov))
  colnames(Z_network) <- make.unique(colnames(Z_network))
  
  W_mean <- cbind(X_cov, Z_network)
  colnames(W_mean) <- make.unique(c(colnames(X_cov), colnames(Z_network)))
  p_mean <- ncol(W_mean)
  mean_groups <- list(
    cov = seq_len(p_cov),
    network = p_cov + seq_len(q_network)
  )
  
  D2_cov <- sqdist(X_cov)
  
  if (is.null(network_geometry)) {
    network_geometry <- prepare_network_components_for_blocks(
      A_network,
      network_name = "network"
    )
  } else {
    if (!is.list(network_geometry) || is.null(network_geometry$A) ||
        !all(dim(network_geometry$A) == c(n, n))) {
      stop("The supplied precomputed network geometry is invalid.")
    }
  }
  
  hp_cov <- hyperpriors$expert0_cov
  hp_network <- hyperpriors$expert1_network
  if (is.null(hp_cov) || is.null(hp_network)) {
    stop("hyperpriors must contain expert0_cov and expert1_network.")
  }
  
  prior_sd_blocks <- list(
    cov = pmax(c(hp_cov$log_ell_sd, hp_cov$log_sig_sd), 1e-8),
    network = pmax(c(
      hp_network$log_ell_sd,
      hp_network$log_sig_sd
    ), 1e-8)
  )
  prior_sd_theta <- unlist(prior_sd_blocks, use.names = FALSE)
  names(prior_sd_theta) <- theta_names
  
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
  
  normalize_proposal_cor_list <- function(x) {
    if (is.null(x)) {
      return(list(cov = diag(2), network = diag(2)))
    }
    if (!is.list(x) || !all(group_names %in% names(x))) {
      stop("proposal_cor_init must contain cov and network.")
    }
    list(
      cov = sanitize_cor_2x2(x$cov, max_abs = cor_max_abs),
      network = sanitize_cor_2x2(x$network, max_abs = cor_max_abs)
    )
  }
  
  pi_spike <- min(max(as.numeric(pi_spike), 0), 1)
  
  K_cov_from_theta <- function(theta) {
    rbf_kernel_from_D2(
      D2_cov,
      sigma_f = exp(theta[2L]),
      l = exp(theta[1L])
    )
  }
  
  K_network_from_theta <- function(theta) {
    network_kernel_from_geom(
      network_geometry,
      ell = exp(theta[1L]),
      sigma = exp(theta[2L])
    )
  }
  
  theta_ok_cov <- function(theta) {
    all(is.finite(theta)) &&
      theta[1L] >= logl_cov_bounds[1L] &&
      theta[1L] <= logl_cov_bounds[2L] &&
      theta[2L] >= min_log_sigma &&
      theta[2L] <= max_log_sigma_cov
  }
  
  theta_ok_network <- function(theta) {
    if (!all(is.finite(theta))) return(FALSE)
    if (theta[1L] < logl_network_bounds[1L] ||
        theta[1L] > logl_network_bounds[2L]) return(FALSE)
    if (theta[2L] < min_log_sigma ||
        theta[2L] > max_log_sigma_network) return(FALSE)
    
    K <- K_network_from_theta(theta)
    mean_marginal_sd <- mean(sqrt(pmax(diag(K), 0)))
    is.finite(mean_marginal_sd) &&
      mean_marginal_sd <= max_mean_sd_network
  }
  
  spike_slab_logsigma <- function(log_sigma, hp) {
    slab_sd <- max(as.numeric(hp$log_sig_sd), 1e-8)
    log_slab <- dnorm(
      log_sigma,
      mean = hp$log_sig_mean,
      sd = slab_sd,
      log = TRUE
    )
    if (pi_spike <= 0) return(log_slab)
    
    spike_sd <- max(as.numeric(spike_log_sig_sd), 1e-8)
    log_spike <- dnorm(
      log_sigma,
      mean = spike_log_sig_mean,
      sd = spike_sd,
      log = TRUE
    )
    if (pi_spike >= 1) return(log_spike)
    
    aa <- log(pi_spike) + log_spike
    bb <- log(1 - pi_spike) + log_slab
    mm <- max(aa, bb)
    mm + log(exp(aa - mm) + exp(bb - mm))
  }
  
  logprior_cov <- function(theta) {
    if (!theta_ok_cov(theta)) return(-Inf)
    dnorm(
      theta[1L],
      mean = hp_cov$log_ell_mean,
      sd = max(hp_cov$log_ell_sd, 1e-8),
      log = TRUE
    ) +
      spike_slab_logsigma(theta[2L], hp_cov)
  }
  
  logprior_network <- function(theta) {
    if (!theta_ok_network(theta)) return(-Inf)
    dnorm(
      theta[1L],
      mean = hp_network$log_ell_mean,
      sd = max(hp_network$log_ell_sd, 1e-8),
      log = TRUE
    ) +
      spike_slab_logsigma(theta[2L], hp_network)
  }
  
  logprior_all <- function(theta_vector) {
    theta_vector <- as.numeric(theta_vector)
    if (length(theta_vector) != 4L) return(-Inf)
    value <- logprior_cov(theta_vector[1:2]) +
      logprior_network(theta_vector[3:4])
    if (!is.finite(value)) return(-Inf)
    value
  }
  
  make_initial_theta_block <- function(hp, kind, override = NULL) {
    theta <- if (is.null(override)) {
      c(hp$log_ell_mean, hp$log_sig_mean)
    } else {
      as.numeric(override)
    }
    if (length(theta) != 2L || any(!is.finite(theta))) {
      stop("Each theta block must contain log_ell and log_sigma.")
    }
    
    if (identical(kind, "cov")) {
      theta[1L] <- min(max(theta[1L], logl_cov_bounds[1L]),
                       logl_cov_bounds[2L])
      theta[2L] <- min(max(theta[2L], min_log_sigma),
                       max_log_sigma_cov)
      if (!theta_ok_cov(theta)) stop("Initial covariate theta is invalid.")
      return(theta)
    }
    
    theta[1L] <- min(max(theta[1L], logl_network_bounds[1L]),
                     logl_network_bounds[2L])
    theta[2L] <- min(max(theta[2L], min_log_sigma),
                     max_log_sigma_network)
    while (!theta_ok_network(theta) && theta[2L] > min_log_sigma) {
      theta[2L] <- max(theta[2L] - log(1.25), min_log_sigma)
    }
    if (!theta_ok_network(theta)) {
      stop("Initial network theta is invalid.")
    }
    theta
  }
  
  Csum_svd <- function(theta_vector) {
    K_cov <- K_cov_from_theta(theta_vector[1:2])
    K_network <- K_network_from_theta(theta_vector[3:4])
    Csum <- symmetrize(K_cov + K_network) + diag(jitter, n)
    svd_psd_fallback(Csum)
  }
  
  Z_mean <- cbind(beta0 = 1, W_mean)
  colnames(Z_mean) <- c(
    "beta0",
    paste0("gamma_", colnames(W_mean))
  )
  p_coef <- ncol(Z_mean)
  prior_mean_coef <- c(m_beta, rep(0, p_mean))
  names(prior_mean_coef) <- colnames(Z_mean)
  
  if (is.null(init_state)) {
    theta_vector <- c(
      make_initial_theta_block(hp_cov, "cov"),
      make_initial_theta_block(hp_network, "network")
    )
    names(theta_vector) <- theta_names
    
    coef_vec <- prior_mean_coef
    tau2_group <- tau2_default
    aux_tau <- setNames(numeric(length(group_names)), group_names)
    for (group_name in group_names) {
      aux_tau[group_name] <- rInvGamma(
        1L,
        shape = 1 / 2,
        scale = 1 / (scale_tau[group_name]^2)
      )
    }
    f_init_supplied <- NULL
  } else {
    theta_input <- as.numeric(init_state$theta_log)
    if (length(theta_input) != 4L || any(!is.finite(theta_input))) {
      stop("The initial theta state is invalid.")
    }
    theta_vector <- c(
      make_initial_theta_block(hp_cov, "cov", theta_input[1:2]),
      make_initial_theta_block(
        hp_network,
        "network",
        theta_input[3:4]
      )
    )
    names(theta_vector) <- theta_names
    
    coef_vec <- as.numeric(init_state$beta)
    tau2_group <- normalize_named_positive(
      init_state$tau2,
      group_names,
      "init_state$tau2"
    )
    aux_tau <- normalize_named_positive(
      init_state$aux_tau,
      group_names,
      "init_state$aux_tau"
    )
    f_init_supplied <- init_state$f
  }
  
  if (length(coef_vec) != p_coef || any(!is.finite(coef_vec))) {
    stop("The initial beta/gamma state has the wrong length or invalid values.")
  }
  names(coef_vec) <- colnames(Z_mean)
  
  mu_vec <- as.vector(Z_mean %*% coef_vec)
  if (is.null(f_init_supplied)) {
    f <- mu_vec
  } else {
    f <- as.numeric(f_init_supplied)
    if (length(f) != n || any(!is.finite(f))) {
      stop("The supplied initial f state is invalid.")
    }
  }
  
  svd_sum <- Csum_svd(theta_vector)
  U_sum <- svd_sum$U
  gam_sum <- pmax(svd_sum$d, 0)
  
  log_likelihood <- bernoulli_loglik_vec(f, Y)
  gradient_f <- bernoulli_gradloglik_vec(f, Y)
  
  kappa_cov <- normalize_scalar(kappa_cov_init, "kappa_cov_init")
  kappa_network <- normalize_scalar(
    kappa_network_init,
    "kappa_network_init"
  )
  proposal_cor <- normalize_proposal_cor_list(proposal_cor_init)
  joint_theta_scale <- normalize_scalar(
    joint_theta_scale_init,
    "joint_theta_scale_init"
  )
  delta_state <- normalize_scalar(delta_init, "delta_init")
  
  prior_var_from_tau <- function(tau2_current) {
    output <- numeric(p_mean)
    output[mean_groups$cov] <- tau2_current["cov"]
    output[mean_groups$network] <- tau2_current["network"]
    pmax(output, .Machine$double.eps)
  }
  
  sample_beta_gamma_gibbs <- function(f, U, gam, tau2_current) {
    prior_var <- c(s_beta^2, prior_var_from_tau(tau2_current))
    eig <- pmax(as.numeric(gam), .Machine$double.eps)
    
    UtZ <- crossprod(U, Z_mean)
    KinvZ <- U %*% sweep(UtZ, 1L, eig, "/")
    Utf <- as.vector(crossprod(U, f))
    Kinvf <- as.vector(U %*% (Utf / eig))
    
    precision <- symmetrize(
      crossprod(Z_mean, KinvZ) + diag(1 / prior_var, p_coef)
    )
    information <- as.vector(crossprod(Z_mean, Kinvf)) +
      prior_mean_coef / prior_var
    
    R <- chol_safe(precision, jitter = 1e-10, max_tries = 8L)$R
    mean_coef <- as.vector(
      backsolve(R, forwardsolve(t(R), information))
    )
    mean_coef + as.vector(backsolve(R, rnorm(p_coef)))
  }
  
  update_coef_tau <- function(
    f,
    U,
    gam,
    coef_current,
    tau2_current,
    aux_current
  ) {
    coef_new <- sample_beta_gamma_gibbs(
      f,
      U,
      gam,
      tau2_current
    )
    names(coef_new) <- colnames(Z_mean)
    
    if (isTRUE(sample_tau_groups)) {
      gamma_new <- coef_new[-1L]
      for (group_name in group_names) {
        idx <- mean_groups[[group_name]]
        tau2_current[group_name] <- rInvGamma(
          1L,
          shape = (nu_tau + length(idx)) / 2,
          scale = (nu_tau / aux_current[group_name]) +
            sum(gamma_new[idx]^2) / 2
        )
        aux_current[group_name] <- rInvGamma(
          1L,
          shape = (nu_tau + 1) / 2,
          scale = (1 / scale_tau[group_name]^2) +
            (nu_tau / tau2_current[group_name])
        )
      }
    }
    
    list(
      coef = coef_new,
      tau2 = tau2_current,
      aux_tau = aux_current,
      mu = as.vector(Z_mean %*% coef_new)
    )
  }
  
  adapt_log_step <- function(
    step,
    acceptance,
    target,
    iteration,
    rate,
    t0,
    power,
    minimum,
    maximum
  ) {
    gain <- rate / ((iteration + t0)^power)
    updated <- step * exp(gain * (acceptance - target))
    max(minimum, min(maximum, updated))
  }
  
  agrad_f_step <- function(
    f,
    mu,
    log_likelihood,
    gradient_f,
    U,
    gam,
    delta
  ) {
    nn <- length(f)
    residual <- f - mu
    z <- residual +
      (delta / 2) * gradient_f +
      sqrt(delta / 2) * rnorm(nn)
    
    sqrt_lambda_tilde <- sqrt(
      (pmax(gam, 0) * delta) /
        pmax(delta + 2 * pmax(gam, 0), .Machine$double.eps)
    )
    sqrt_lambda_tilde[!is.finite(sqrt_lambda_tilde)] <- 0
    
    temp1 <- crossprod(U, (2 / delta) * z)
    temp2 <- sqrt_lambda_tilde * temp1 + rnorm(nn)
    temp3 <- sqrt_lambda_tilde * temp2
    residual_proposed <- as.vector(U %*% temp3)
    f_proposed <- mu + residual_proposed
    
    log_likelihood_proposed <- bernoulli_loglik_vec(f_proposed, Y)
    gradient_proposed <- bernoulli_gradloglik_vec(f_proposed, Y)
    
    gzy <- as.numeric(crossprod(
      z - residual_proposed - (delta / 4) * gradient_proposed,
      gradient_proposed
    ))
    gzx <- as.numeric(crossprod(
      z - residual - (delta / 4) * gradient_f,
      gradient_f
    ))
    
    log_alpha <-
      (log_likelihood_proposed - log_likelihood) +
      (gzy - gzx)
    accepted <- log(runif(1)) < min(0, log_alpha)
    
    if (accepted) {
      list(
        f = f_proposed,
        log_likelihood = log_likelihood_proposed,
        gradient_f = gradient_proposed,
        accepted = TRUE
      )
    } else {
      list(
        f = f,
        log_likelihood = log_likelihood,
        gradient_f = gradient_f,
        accepted = FALSE
      )
    }
  }
  
  theta_block_step <- function(
    block,
    theta_current,
    f_fixed,
    mu_fixed,
    U_current,
    gam_current,
    kappa_scalar,
    correlation_block,
    prior_sd_block
  ) {
    log_prior_current <- logprior_all(theta_current)
    log_normal_current <- logN_mu_eig(
      f_fixed,
      mu_fixed,
      U_current,
      gam_current
    )
    
    proposal_covariance <- proposal_cov_scaled_2x2(
      kappa_scalar = kappa_scalar,
      prior_sd = prior_sd_block,
      cor_mat = correlation_block
    )
    proposal_increment <- rmvnorm0_chol(proposal_covariance)
    
    theta_proposed <- theta_current
    if (identical(block, "cov")) {
      theta_proposed[1:2] <- theta_proposed[1:2] + proposal_increment
    } else if (identical(block, "network")) {
      theta_proposed[3:4] <- theta_proposed[3:4] + proposal_increment
    } else {
      stop("Unknown theta block: ", block)
    }
    
    log_prior_proposed <- logprior_all(theta_proposed)
    if (!is.finite(log_prior_proposed)) {
      return(list(
        theta = theta_current,
        U = U_current,
        gam = gam_current,
        accepted = FALSE
      ))
    }
    
    svd_proposed <- Csum_svd(theta_proposed)
    U_proposed <- svd_proposed$U
    gam_proposed <- pmax(svd_proposed$d, 0)
    log_normal_proposed <- logN_mu_eig(
      f_fixed,
      mu_fixed,
      U_proposed,
      gam_proposed
    )
    
    log_alpha <-
      (log_normal_proposed - log_normal_current) +
      (log_prior_proposed - log_prior_current)
    accepted <- log(runif(1)) < min(0, log_alpha)
    
    if (accepted) {
      list(
        theta = theta_proposed,
        U = U_proposed,
        gam = gam_proposed,
        accepted = TRUE
      )
    } else {
      list(
        theta = theta_current,
        U = U_current,
        gam = gam_current,
        accepted = FALSE
      )
    }
  }
  
  joint_step <- function(
    f,
    mu,
    log_likelihood,
    gradient_f,
    theta_current,
    U_current,
    gam_current,
    delta,
    theta_covariance
  ) {
    nn <- length(f)
    residual <- f - mu
    z <- residual +
      (delta / 2) * gradient_f +
      sqrt(delta / 2) * rnorm(nn)
    
    theta_proposed <- theta_current + rmvnorm0_chol(theta_covariance)
    log_prior_current <- logprior_all(theta_current)
    log_prior_proposed <- logprior_all(theta_proposed)
    
    if (!is.finite(log_prior_proposed)) {
      return(list(
        f = f,
        log_likelihood = log_likelihood,
        gradient_f = gradient_f,
        theta = theta_current,
        U = U_current,
        gam = gam_current,
        accepted = FALSE
      ))
    }
    
    svd_proposed <- Csum_svd(theta_proposed)
    U_proposed <- svd_proposed$U
    gam_proposed <- pmax(svd_proposed$d, 0)
    
    sqrt_lambda_tilde <- sqrt(
      (pmax(gam_proposed, 0) * delta) /
        pmax(delta + 2 * pmax(gam_proposed, 0), .Machine$double.eps)
    )
    sqrt_lambda_tilde[!is.finite(sqrt_lambda_tilde)] <- 0
    
    temp1 <- crossprod(U_proposed, (2 / delta) * z)
    temp2 <- sqrt_lambda_tilde * temp1 + rnorm(nn)
    temp3 <- sqrt_lambda_tilde * temp2
    residual_proposed <- as.vector(U_proposed %*% temp3)
    f_proposed <- mu + residual_proposed
    
    log_likelihood_proposed <- bernoulli_loglik_vec(f_proposed, Y)
    gradient_proposed <- bernoulli_gradloglik_vec(f_proposed, Y)
    
    gzy <- as.numeric(crossprod(
      z - residual_proposed - (delta / 4) * gradient_proposed,
      gradient_proposed
    ))
    gzx <- as.numeric(crossprod(
      z - residual - (delta / 4) * gradient_f,
      gradient_f
    ))
    
    logZ_current <- logN0_eig(
      z,
      U_current,
      gam_current,
      nugget = delta / 2,
      jitter = jitter
    )
    logZ_proposed <- logN0_eig(
      z,
      U_proposed,
      gam_proposed,
      nugget = delta / 2,
      jitter = jitter
    )
    
    log_alpha <-
      (log_likelihood_proposed - log_likelihood) +
      (gzy - gzx) +
      (logZ_proposed - logZ_current) +
      (log_prior_proposed - log_prior_current)
    accepted <- log(runif(1)) < min(0, log_alpha)
    
    if (accepted) {
      list(
        f = f_proposed,
        log_likelihood = log_likelihood_proposed,
        gradient_f = gradient_proposed,
        theta = theta_proposed,
        U = U_proposed,
        gam = gam_proposed,
        accepted = TRUE
      )
    } else {
      list(
        f = f,
        log_likelihood = log_likelihood,
        gradient_f = gradient_f,
        theta = theta_current,
        U = U_current,
        gam = gam_current,
        accepted = FALSE
      )
    }
  }
  
  ## Pilot tail storage. Only cross-data-set-comparable coefficients are kept:
  ## beta0 and the ten covariate coefficients. Network coefficient coordinates
  ## are not pooled across independently generated networks.
  tail_theta <- if (state_tail_length > 0L) {
    matrix(
      NA_real_,
      nrow = state_tail_length,
      ncol = 4L,
      dimnames = list(NULL, theta_names)
    )
  } else {
    NULL
  }
  tail_beta_fixed <- if (state_tail_length > 0L) {
    matrix(
      NA_real_,
      nrow = state_tail_length,
      ncol = 11L,
      dimnames = list(NULL, colnames(Z_mean)[1:11])
    )
  } else {
    NULL
  }
  tail_tau2 <- if (state_tail_length > 0L) {
    matrix(
      NA_real_,
      nrow = state_tail_length,
      ncol = 2L,
      dimnames = list(NULL, group_names)
    )
  } else {
    NULL
  }
  tail_aux_tau <- if (state_tail_length > 0L) {
    matrix(
      NA_real_,
      nrow = state_tail_length,
      ncol = 2L,
      dimnames = list(NULL, group_names)
    )
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
      tail_theta[tail_index, ] <<- theta_vector
      tail_beta_fixed[tail_index, ] <<- coef_vec[1:11]
      tail_tau2[tail_index, ] <<- tau2_group[group_names]
      tail_aux_tau[tail_index, ] <<- aux_tau[group_names]
    }
    invisible(NULL)
  }
  
  start_time <- Sys.time()
  if (verbose) {
    cat(sprintf(
      "\n[%s] %s exact two-expert sampler started.\n",
      timestamp_now(),
      phase_label
    ))
    cat("Observations            :", n, "\n")
    cat("Covariate columns       :", p_cov, "\n")
    cat("Network mean columns    :", q_network, "\n")
    cat("Blocked warm-up         :", blocked_iters, "\n")
    cat("Joint tuning            :", joint_tune_iters, "\n")
    cat("Posterior sampling      :", n_sample, "\n")
    cat("Thinning                :", thin, "\n")
    cat("Retained draws          :", n_saved_expected, "\n")
    cat("Saved pilot tail        :", state_tail_length, "\n")
  }
  
  ## ------------------------------------------------------------------------
  ## Phase 1: correlated blocked warm-up
  ## ------------------------------------------------------------------------
  
  acceptance_blocked_sum <- c(f = 0L, cov = 0L, network = 0L)
  adaptation_window_length <- 0L
  adaptation_window_acceptance <- c(f = 0L, cov = 0L, network = 0L)
  adaptation_window_id <- 0L
  
  theta_path_blocked <- matrix(
    NA_real_,
    nrow = max(blocked_iters, 1L),
    ncol = 4L,
    dimnames = list(NULL, theta_names)
  )
  cor_adapt_trace <- data.frame(
    iteration = integer(0),
    rho_cov = numeric(0),
    rho_network = numeric(0)
  )
  
  if (blocked_iters > 0L) {
    for (iteration in seq_len(blocked_iters)) {
      f_update <- agrad_f_step(
        f,
        mu_vec,
        log_likelihood,
        gradient_f,
        U_sum,
        gam_sum,
        delta_state
      )
      f <- f_update$f
      log_likelihood <- f_update$log_likelihood
      gradient_f <- f_update$gradient_f
      accepted_f <- as.integer(f_update$accepted)
      
      theta_update <- theta_block_step(
        block = "cov",
        theta_current = theta_vector,
        f_fixed = f,
        mu_fixed = mu_vec,
        U_current = U_sum,
        gam_current = gam_sum,
        kappa_scalar = kappa_cov,
        correlation_block = proposal_cor$cov,
        prior_sd_block = prior_sd_blocks$cov
      )
      theta_vector <- theta_update$theta
      U_sum <- theta_update$U
      gam_sum <- theta_update$gam
      accepted_cov <- as.integer(theta_update$accepted)
      
      theta_update <- theta_block_step(
        block = "network",
        theta_current = theta_vector,
        f_fixed = f,
        mu_fixed = mu_vec,
        U_current = U_sum,
        gam_current = gam_sum,
        kappa_scalar = kappa_network,
        correlation_block = proposal_cor$network,
        prior_sd_block = prior_sd_blocks$network
      )
      theta_vector <- theta_update$theta
      U_sum <- theta_update$U
      gam_sum <- theta_update$gam
      accepted_network <- as.integer(theta_update$accepted)
      
      remaining_update <- update_coef_tau(
        f,
        U_sum,
        gam_sum,
        coef_vec,
        tau2_group,
        aux_tau
      )
      coef_vec <- remaining_update$coef
      tau2_group <- remaining_update$tau2
      aux_tau <- remaining_update$aux_tau
      mu_vec <- remaining_update$mu
      
      acceptance_now <- c(
        f = accepted_f,
        cov = accepted_cov,
        network = accepted_network
      )
      acceptance_blocked_sum <- acceptance_blocked_sum + acceptance_now
      theta_path_blocked[iteration, ] <- theta_vector
      
      if (adapt && adapt_block_cor &&
          iteration >= cor_adapt_start &&
          iteration %% cor_adapt_interval == 0L) {
        idx_cor <- seq.int(
          max(1L, iteration - cor_adapt_window + 1L),
          iteration
        )
        proposal_cor$cov <- estimate_cor_2x2(
          theta_path_blocked[idx_cor, 1:2, drop = FALSE],
          fallback = proposal_cor$cov,
          shrinkage = cor_shrinkage,
          max_abs = cor_max_abs,
          min_n = min(200L, length(idx_cor))
        )
        proposal_cor$network <- estimate_cor_2x2(
          theta_path_blocked[idx_cor, 3:4, drop = FALSE],
          fallback = proposal_cor$network,
          shrinkage = cor_shrinkage,
          max_abs = cor_max_abs,
          min_n = min(200L, length(idx_cor))
        )
        cor_adapt_trace <- rbind(
          cor_adapt_trace,
          data.frame(
            iteration = iteration,
            rho_cov = proposal_cor$cov[1L, 2L],
            rho_network = proposal_cor$network[1L, 2L]
          )
        )
      }
      
      if (adapt && iteration >= adapt_start) {
        adaptation_window_length <- adaptation_window_length + 1L
        adaptation_window_acceptance <-
          adaptation_window_acceptance + acceptance_now
        
        if (adaptation_window_length == adapt_interval) {
          adaptation_window_id <- adaptation_window_id + 1L
          acceptance_rate <-
            adaptation_window_acceptance / adaptation_window_length
          
          delta_state <- adapt_log_step(
            delta_state,
            acceptance_rate["f"],
            target_f,
            adaptation_window_id,
            adapt_rate_delta,
            t0_adapt,
            power_adapt,
            delta_min,
            delta_max
          )
          kappa_cov <- adapt_log_step(
            kappa_cov,
            acceptance_rate["cov"],
            target_theta,
            adaptation_window_id,
            adapt_rate_theta,
            t0_adapt,
            power_adapt,
            kappa_min,
            kappa_max
          )
          kappa_network <- adapt_log_step(
            kappa_network,
            acceptance_rate["network"],
            target_theta,
            adaptation_window_id,
            adapt_rate_theta,
            t0_adapt,
            power_adapt,
            kappa_min,
            kappa_max
          )
          
          adaptation_window_length <- 0L
          adaptation_window_acceptance[] <- 0L
        }
      }
      
      store_tail_state(iteration)
      
      if (verbose && burnin_progress_every > 0L &&
          iteration %% burnin_progress_every == 0L) {
        progress_message(
          paste0(phase_label, " blocked warm-up"),
          iteration,
          blocked_iters
        )
      }
    }
  }
  
  acceptance_blocked <-
    acceptance_blocked_sum / max(1L, blocked_iters)
  proposal_cor_final <- list(
    cov = sanitize_cor_2x2(proposal_cor$cov, max_abs = cor_max_abs),
    network = sanitize_cor_2x2(
      proposal_cor$network,
      max_abs = cor_max_abs
    )
  )
  
  ## ------------------------------------------------------------------------
  ## Construct the full 4 x 4 joint proposal geometry
  ## ------------------------------------------------------------------------
  
  fallback_correlation <- diag(4L)
  fallback_correlation[1:2, 1:2] <- proposal_cor_final$cov
  fallback_correlation[3:4, 3:4] <- proposal_cor_final$network
  dimnames(fallback_correlation) <- list(theta_names, theta_names)
  
  theta_tail_blocked <- if (blocked_iters > 0L) {
    idx_start <- max(1L, blocked_iters - emp_cov_tail + 1L)
    theta_path_blocked[idx_start:blocked_iters, , drop = FALSE]
  } else {
    matrix(theta_vector, nrow = 1L, dimnames = list(NULL, theta_names))
  }
  
  if (!is.null(joint_cov_base_init)) {
    Sigma_theta_base <- sanitize_covariance_matrix(
      as.matrix(joint_cov_base_init) + cov_jitter * diag(4L)
    )
    if (!all(dim(Sigma_theta_base) == c(4L, 4L))) {
      stop("joint_cov_base_init must be 4 x 4.")
    }
    joint_correlation_base <- sanitize_correlation_matrix(
      safe_cor_from_cov(Sigma_theta_base),
      max_abs = joint_cor_max_abs,
      shrinkage = 0
    )
    joint_geometry_source <- "supplied_global_pilot_covariance"
  } else {
    joint_correlation_base <- estimate_correlation_matrix(
      theta_tail_blocked,
      fallback = fallback_correlation,
      shrinkage = joint_cor_shrinkage,
      max_abs = joint_cor_max_abs,
      min_n = min(200L, nrow(theta_tail_blocked))
    )
    
    block_scale <- c(
      rep(sqrt(kappa_cov), 2L),
      rep(sqrt(kappa_network), 2L)
    )
    base_sd <- prior_sd_theta * block_scale
    D_theta <- diag(base_sd, 4L)
    Sigma_theta_base <- D_theta %*%
      joint_correlation_base %*%
      D_theta
    Sigma_theta_base <- sanitize_covariance_matrix(
      Sigma_theta_base + cov_jitter * diag(4L)
    )
    joint_geometry_source <-
      "blocked_path_prior_scaled_full_correlation"
  }
  
  dimnames(joint_correlation_base) <- list(theta_names, theta_names)
  dimnames(Sigma_theta_base) <- list(theta_names, theta_names)
  
  delta_joint <- delta_state
  
  if (verbose) {
    cat(sprintf(
      "[%s] %s blocked phase complete.\n",
      timestamp_now(),
      phase_label
    ))
    cat("  delta_joint          :", delta_joint, "\n")
    cat("  kappa_cov            :", kappa_cov, "\n")
    cat("  kappa_network        :", kappa_network, "\n")
    cat("  joint geometry source:", joint_geometry_source, "\n")
  }
  
  ## ------------------------------------------------------------------------
  ## Phase 2: joint tuning
  ## ------------------------------------------------------------------------
  
  acceptance_joint_tune_sum <- 0L
  acceptance_f_refresh_tune_sum <- 0L
  joint_window_length <- 0L
  joint_window_acceptance <- 0L
  joint_window_id <- 0L
  
  if (joint_tune_iters > 0L) {
    for (iteration in seq_len(joint_tune_iters)) {
      Sigma_theta_effective <-
        (joint_theta_scale^2) * Sigma_theta_base
      
      joint_update <- joint_step(
        f,
        mu_vec,
        log_likelihood,
        gradient_f,
        theta_vector,
        U_sum,
        gam_sum,
        delta_joint,
        Sigma_theta_effective
      )
      f <- joint_update$f
      log_likelihood <- joint_update$log_likelihood
      gradient_f <- joint_update$gradient_f
      theta_vector <- joint_update$theta
      U_sum <- joint_update$U
      gam_sum <- joint_update$gam
      accepted_joint <- as.integer(joint_update$accepted)
      acceptance_joint_tune_sum <-
        acceptance_joint_tune_sum + accepted_joint
      
      if (isTRUE(f_refresh_after_joint)) {
        f_update <- agrad_f_step(
          f,
          mu_vec,
          log_likelihood,
          gradient_f,
          U_sum,
          gam_sum,
          delta_joint
        )
        f <- f_update$f
        log_likelihood <- f_update$log_likelihood
        gradient_f <- f_update$gradient_f
        acceptance_f_refresh_tune_sum <-
          acceptance_f_refresh_tune_sum +
          as.integer(f_update$accepted)
      }
      
      remaining_update <- update_coef_tau(
        f,
        U_sum,
        gam_sum,
        coef_vec,
        tau2_group,
        aux_tau
      )
      coef_vec <- remaining_update$coef
      tau2_group <- remaining_update$tau2
      aux_tau <- remaining_update$aux_tau
      mu_vec <- remaining_update$mu
      
      if (adapt && iteration >= adapt_start) {
        joint_window_length <- joint_window_length + 1L
        joint_window_acceptance <-
          joint_window_acceptance + accepted_joint
        
        if (joint_window_length == adapt_interval) {
          joint_window_id <- joint_window_id + 1L
          joint_rate <- joint_window_acceptance / joint_window_length
          joint_theta_scale <- adapt_log_step(
            joint_theta_scale,
            joint_rate,
            target_joint,
            joint_window_id,
            adapt_rate_joint_theta,
            t0_adapt,
            power_adapt,
            joint_theta_scale_min,
            joint_theta_scale_max
          )
          joint_window_length <- 0L
          joint_window_acceptance <- 0L
        }
      }
      
      store_tail_state(blocked_iters + iteration)
      
      if (verbose && burnin_progress_every > 0L &&
          iteration %% burnin_progress_every == 0L) {
        progress_message(
          paste0(phase_label, " joint tuning"),
          iteration,
          joint_tune_iters
        )
      }
    }
  }
  
  acceptance_joint_tune <-
    acceptance_joint_tune_sum / max(1L, joint_tune_iters)
  acceptance_f_refresh_tune <-
    acceptance_f_refresh_tune_sum / max(1L, joint_tune_iters)
  
  delta_final <- delta_joint
  joint_theta_scale_final <- joint_theta_scale
  Sigma_theta_final <- sanitize_covariance_matrix(
    (joint_theta_scale_final^2) * Sigma_theta_base
  )
  dimnames(Sigma_theta_final) <- list(theta_names, theta_names)
  
  block_covariance_final <- list(
    cov = sanitize_covariance_matrix(
      proposal_cov_scaled_2x2(
        kappa_cov,
        prior_sd_blocks$cov,
        proposal_cor_final$cov,
        jitter = cov_jitter
      )
    ),
    network = sanitize_covariance_matrix(
      proposal_cov_scaled_2x2(
        kappa_network,
        prior_sd_blocks$network,
        proposal_cor_final$network,
        jitter = cov_jitter
      )
    )
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
  
  ## ------------------------------------------------------------------------
  ## Phase 3: frozen posterior sampling
  ## ------------------------------------------------------------------------
  
  score_sum <- numeric(n)
  n_saved <- 0L
  
  predictive_label_probability_draws <- if (n_saved_expected > 0L) {
    matrix(
      NA_real_,
      nrow = n_saved_expected,
      ncol = n,
      dimnames = list(
        NULL,
        paste0("node_", seq_len(n))
      )
    )
  } else {
    NULL
  }
  
  max_probability_draws <- numeric(n_saved_expected)
  q99_probability_draws <- numeric(n_saved_expected)
  top5_mean_probability_draws <- numeric(n_saved_expected)
  beta0_draws <- numeric(n_saved_expected)
  tau_draws <- if (n_saved_expected > 0L) {
    matrix(
      NA_real_,
      nrow = n_saved_expected,
      ncol = 2L,
      dimnames = list(NULL, c("tau_cov", "tau_network"))
    )
  } else {
    matrix(numeric(0), nrow = 0L, ncol = 2L)
  }
  theta_draws <- if (n_saved_expected > 0L) {
    matrix(
      NA_real_,
      nrow = n_saved_expected,
      ncol = 4L,
      dimnames = list(
        NULL,
        c("ell_cov", "sigma_cov", "ell_network", "sigma_network")
      )
    )
  } else {
    matrix(numeric(0), nrow = 0L, ncol = 4L)
  }
  
  acceptance_joint_sampling_sum <- 0L
  acceptance_f_refresh_sampling_sum <- 0L
  
  if (n_sample > 0L) {
    for (iteration in seq_len(n_sample)) {
      joint_update <- joint_step(
        f,
        mu_vec,
        log_likelihood,
        gradient_f,
        theta_vector,
        U_sum,
        gam_sum,
        delta_final,
        Sigma_theta_final
      )
      f <- joint_update$f
      log_likelihood <- joint_update$log_likelihood
      gradient_f <- joint_update$gradient_f
      theta_vector <- joint_update$theta
      U_sum <- joint_update$U
      gam_sum <- joint_update$gam
      acceptance_joint_sampling_sum <-
        acceptance_joint_sampling_sum +
        as.integer(joint_update$accepted)
      
      if (isTRUE(f_refresh_after_joint)) {
        f_update <- agrad_f_step(
          f,
          mu_vec,
          log_likelihood,
          gradient_f,
          U_sum,
          gam_sum,
          delta_final
        )
        f <- f_update$f
        log_likelihood <- f_update$log_likelihood
        gradient_f <- f_update$gradient_f
        acceptance_f_refresh_sampling_sum <-
          acceptance_f_refresh_sampling_sum +
          as.integer(f_update$accepted)
      }
      
      remaining_update <- update_coef_tau(
        f,
        U_sum,
        gam_sum,
        coef_vec,
        tau2_group,
        aux_tau
      )
      coef_vec <- remaining_update$coef
      tau2_group <- remaining_update$tau2
      aux_tau <- remaining_update$aux_tau
      mu_vec <- remaining_update$mu
      
      if (iteration %% thin == 0L) {
        n_saved <- n_saved + 1L
        p_current <- plogis(f)
        score_sum <- score_sum + p_current
        
        predictive_label_probability_draws[
          n_saved,
          seq_len(n)
        ] <- p_current
        
        max_probability_draws[n_saved] <- max(p_current)
        q99_probability_draws[n_saved] <- unname(
          quantile(p_current, 0.99, names = FALSE)
        )
        top5_mean_probability_draws[n_saved] <- mean(
          sort(p_current, decreasing = TRUE)[1:5]
        )
        beta0_draws[n_saved] <- coef_vec[1L]
        tau_draws[n_saved, ] <- sqrt(tau2_group[group_names])
        theta_draws[n_saved, ] <- exp(theta_vector)
      }
      
      if (verbose && sample_progress_every > 0L &&
          iteration %% sample_progress_every == 0L) {
        progress_message(
          paste0(phase_label, " posterior sampling"),
          iteration,
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
  
  predictive_T_draws <- matrix(
    integer(0),
    nrow = 0L,
    ncol = n
  )
  if (n_saved_expected > 0L) {
    if (is.null(predictive_label_seed)) {
      stop(
        "predictive_label_seed is required when posterior draws are retained."
      )
    }
    if (any(!is.finite(predictive_label_probability_draws))) {
      stop("The predictive label probability matrix was not filled completely.")
    }
    predictive_T_draws <- simulate_binary_label_draws(
      predictive_label_probability_draws,
      seed = predictive_label_seed
    )
    rm(predictive_label_probability_draws)
  }
  
  end_time <- Sys.time()
  runtime_sec <- as.numeric(difftime(
    end_time,
    start_time,
    units = "secs"
  ))
  
  acceptance_joint_sampling <-
    acceptance_joint_sampling_sum / max(1L, n_sample)
  acceptance_f_refresh_sampling <-
    acceptance_f_refresh_sampling_sum / max(1L, n_sample)
  
  if (verbose) {
    cat(sprintf(
      "[%s] %s exact two-expert sampler finished.\n",
      timestamp_now(),
      phase_label
    ))
    cat("Runtime                    :", format_duration(runtime_sec), "\n")
    cat("Blocked acceptance f       :", acceptance_blocked["f"], "\n")
    cat("Blocked acceptance cov     :", acceptance_blocked["cov"], "\n")
    cat("Blocked acceptance network :", acceptance_blocked["network"], "\n")
    cat("Joint tune acceptance      :", acceptance_joint_tune, "\n")
    cat("Joint sample acceptance    :", acceptance_joint_sampling, "\n")
    cat("f-refresh sample accept.   :", acceptance_f_refresh_sampling, "\n")
    cat("Final ell/sigma cov        :",
        exp(theta_vector[1L]), "/", exp(theta_vector[2L]), "\n")
    cat("Final ell/sigma network    :",
        exp(theta_vector[3L]), "/", exp(theta_vector[4L]), "\n")
    cat("Final tau values           :\n")
    print(sqrt(tau2_group))
    cat("Final delta                :", delta_final, "\n")
    cat("Final joint theta scale    :", joint_theta_scale_final, "\n")
  }
  
  list(
    state = list(
      theta_log = setNames(as.numeric(theta_vector), theta_names),
      beta = setNames(as.numeric(coef_vec), colnames(Z_mean)),
      tau2 = setNames(as.numeric(tau2_group), group_names),
      aux_tau = setNames(as.numeric(aux_tau), group_names),
      f = as.numeric(f)
    ),
    tail = if (state_tail_length > 0L) {
      list(
        theta_log = tail_theta,
        beta_fixed = tail_beta_fixed,
        tau2 = tail_tau2,
        aux_tau = tail_aux_tau
      )
    } else {
      NULL
    },
    posterior_mean_scores = if (n_saved > 0L) {
      score_sum / n_saved
    } else {
      NULL
    },
    predictive_T_draws = predictive_T_draws,
    max_probability_draws = max_probability_draws,
    q99_probability_draws = q99_probability_draws,
    top5_mean_probability_draws = top5_mean_probability_draws,
    beta0_draws = beta0_draws,
    tau_draws = tau_draws,
    theta_draws = theta_draws,
    n_saved = n_saved,
    network_feature_count = q_network,
    network_component_summary = network_geometry$component_summary,
    proposal = list(
      delta_final = as.numeric(delta_final),
      kappa_cov_final = as.numeric(kappa_cov),
      kappa_network_final = as.numeric(kappa_network),
      proposal_cor_final = proposal_cor_final,
      block_covariance_final = block_covariance_final,
      joint_correlation_base = joint_correlation_base,
      Sigma_theta_base = Sigma_theta_base,
      joint_theta_scale_final = as.numeric(joint_theta_scale_final),
      Sigma_theta_final = Sigma_theta_final,
      joint_geometry_source = joint_geometry_source,
      cor_adapt_trace = cor_adapt_trace
    ),
    acceptance = list(
      blocked_f = unname(acceptance_blocked["f"]),
      blocked_cov = unname(acceptance_blocked["cov"]),
      blocked_network = unname(acceptance_blocked["network"]),
      joint_tune = acceptance_joint_tune,
      f_refresh_tune = acceptance_f_refresh_tune,
      joint_sampling = acceptance_joint_sampling,
      f_refresh_sampling = acceptance_f_refresh_sampling
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
## 11. FEATURE, GEOMETRY, LABEL AND EXACT-PRIOR BUNDLES
## ============================================================================

validate_sar_label_bundle <- function(
    observed_y,
    hidden_positive_idx,
    observed_positive_idx,
    signal_type
) {
  observed_y <- as.integer(observed_y)
  hidden_positive_idx <- sort(unique(as.integer(hidden_positive_idx)))
  observed_positive_idx <- sort(unique(as.integer(observed_positive_idx)))
  signal_type <- as.character(signal_type)

  if (
    length(observed_y) != N_NODES ||
    !all(observed_y %in% c(0L, 1L)) ||
    sum(observed_y) != N_OBSERVED_POS
  ) {
    stop("A feature bundle contains an invalid SAR observed-label vector.")
  }

  if (
    length(hidden_positive_idx) != N_HIDDEN ||
    length(observed_positive_idx) != N_OBSERVED_POS ||
    any(!hidden_positive_idx %in% TRUE_POS_IDX) ||
    any(!observed_positive_idx %in% TRUE_POS_IDX) ||
    length(intersect(hidden_positive_idx, observed_positive_idx)) > 0L ||
    !identical(
      sort(c(hidden_positive_idx, observed_positive_idx)),
      TRUE_POS_IDX
    ) ||
    any(observed_y[hidden_positive_idx] != 0L) ||
    any(observed_y[observed_positive_idx] != 1L)
  ) {
    stop("A feature bundle contains inconsistent SAR positive-label indices.")
  }

  if (length(signal_type) != N_NODES || any(is.na(signal_type))) {
    stop("A feature bundle contains an invalid signal_type vector.")
  }

  list(
    observed_y = observed_y,
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = observed_positive_idx,
    signal_type = signal_type
  )
}

finalize_feature_bundle <- function(
    source_label,
    X_raw,
    X_cov,
    A_network,
    Z_network,
    network_summary,
    modularity_q,
    modularity_selected_eigenvalues,
    seeds,
    observed_y,
    hidden_positive_idx,
    observed_positive_idx,
    signal_type,
    sar_selection_summary,
    sar_selection_diagnostics,
    covariate_profile = NULL,
    X_center = NULL,
    X_scale = NULL,
    calibration_label
) {
  X_raw <- as.matrix(X_raw)
  X_cov <- as.matrix(X_cov)
  A_network <- sanitize_network_adjacency(A_network)
  Z_network <- as.matrix(Z_network)

  labels <- validate_sar_label_bundle(
    observed_y = observed_y,
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = observed_positive_idx,
    signal_type = signal_type
  )

  ## As in the application, every retained modularity
  ## eigenvector is centered and sample-standardized within the current
  ## replication before it enters the structured mean.
  Z_network <- standardize_columns(Z_network)$scaled

  if (
    nrow(X_raw) != N_NODES ||
    ncol(X_raw) != 10L ||
    nrow(X_cov) != N_NODES ||
    ncol(X_cov) != 10L ||
    !all(dim(A_network) == c(N_NODES, N_NODES)) ||
    nrow(Z_network) != N_NODES ||
    ncol(Z_network) < 1L ||
    any(!is.finite(X_raw)) ||
    any(!is.finite(X_cov)) ||
    any(!is.finite(Z_network))
  ) {
    stop("A feature bundle has incompatible dimensions or non-finite values.")
  }

  if (is.null(colnames(X_cov))) {
    colnames(X_cov) <- paste0("X", seq_len(ncol(X_cov)))
  }
  if (is.null(colnames(Z_network))) {
    colnames(Z_network) <- paste0(
      "network_mod",
      seq_len(ncol(Z_network))
    )
  }
  colnames(X_cov) <- make.unique(colnames(X_cov))
  colnames(Z_network) <- make.unique(colnames(Z_network))

  network_geometry <- prepare_network_components_for_blocks(
    A_network,
    network_name = "network"
  )

  prior_calibration <- calibrate_geometry_priors(
    X_cov = X_cov,
    Z_network = Z_network,
    network_geometry = network_geometry,
    calibration_label = calibration_label
  )

  list(
    source = source_label,
    X_raw = X_raw,
    X_cov = X_cov,
    A_network = A_network,
    Z_network = Z_network,
    network_geometry = network_geometry,
    prior_calibration = prior_calibration,
    network_summary = network_summary,
    modularity_q = as.integer(modularity_q),
    modularity_selected_eigenvalues = modularity_selected_eigenvalues,
    seeds = seeds,
    observed_y = labels$observed_y,
    hidden_positive_idx = labels$hidden_positive_idx,
    observed_positive_idx = labels$observed_positive_idx,
    signal_type = labels$signal_type,
    sar_selection_summary = sar_selection_summary,
    sar_selection_diagnostics = sar_selection_diagnostics,
    covariate_profile = covariate_profile,
    X_center = X_center,
    X_scale = X_scale
  )
}

make_generated_feature_bundle <- function(
    seed_row,
    source_label,
    calibration_label
) {
  required_seed_names <- c(
    "covariate_seed",
    "network_seed",
    "modularity_seed",
    "sar_seed"
  )

  if (!all(required_seed_names %in% names(seed_row))) {
    stop("The generated feature bundle is missing required deterministic seeds.")
  }

  covariate_set <- make_covariate_set(
    seed = seed_row$covariate_seed,
    signal_idx = COVARIATE_SIGNAL_IDX
  )

  network_set <- generate_network_set(
    seed = seed_row$network_seed,
    signal_idx = NETWORK_SIGNAL_IDX
  )

  modularity <- modularity_features_single_network(
    A = network_set$A,
    seed = seed_row$modularity_seed,
    eig_tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q,
    k_max = MOD_K_MAX
  )

  sar_design <- build_sar_label_design(
    X_raw = covariate_set$X_raw,
    sar_seed = seed_row$sar_seed
  )

  finalize_feature_bundle(
    source_label = source_label,
    X_raw = covariate_set$X_raw,
    X_cov = covariate_set$X_std,
    A_network = network_set$A,
    Z_network = modularity$Z,
    network_summary = network_set$summary,
    modularity_q = modularity$q,
    modularity_selected_eigenvalues = modularity$selected_eigenvalues,
    seeds = seed_row,
    observed_y = sar_design$observed_y,
    hidden_positive_idx = sar_design$hidden_positive_idx,
    observed_positive_idx = sar_design$observed_positive_idx,
    signal_type = sar_design$signal_type,
    sar_selection_summary = sar_design$selection_summary,
    sar_selection_diagnostics = sar_design$selection_diagnostics,
    covariate_profile = covariate_set$covariate_profile,
    X_center = covariate_set$center,
    X_scale = covariate_set$scale,
    calibration_label = calibration_label
  )
}

load_or_generate_production_bundle <- function(replicate_number) {
  seed_row <- production_seed_table[
    production_seed_table$replicate == replicate_number,
    ,
    drop = FALSE
  ]

  if (nrow(seed_row) != 1L) {
    stop("Production seed lookup failed for replicate ", replicate_number, ".")
  }

  calibration_label <- sprintf("production_%03d", replicate_number)

  geometry_file <- file.path(
    COMPETING_GEOMETRY_DIR,
    sprintf("replicate_%03d_geometry.rds", replicate_number)
  )

  if (isTRUE(USE_SAVED_COMPETING_GEOMETRIES) && file.exists(geometry_file)) {
    geometry <- readRDS(geometry_file)

    required_fields <- c(
      "T",
      "Y",
      "X_raw",
      "X_std",
      "A_network",
      "hidden_positive_idx",
      "observed_positive_idx",
      "signal_type"
    )

    missing_fields <- required_fields[
      !vapply(
        required_fields,
        function(name_now) !is.null(geometry[[name_now]]),
        logical(1)
      )
    ]

    if (length(missing_fields) > 0L) {
      stop(
        "Saved competing SAR geometry is missing: ",
        paste(missing_fields, collapse = ", ")
      )
    }

    if (!identical(as.integer(geometry$T), TRUE_T)) {
      stop("Saved competing SAR geometry has incompatible latent truth T.")
    }

    saved_labels <- validate_sar_label_bundle(
      observed_y = geometry$Y,
      hidden_positive_idx = geometry$hidden_positive_idx,
      observed_positive_idx = geometry$observed_positive_idx,
      signal_type = geometry$signal_type
    )

    if (!is.null(geometry$seeds)) {
      seed_columns <- c(
        "covariate_seed",
        "network_seed",
        "modularity_seed",
        "sar_seed"
      )

      if (all(seed_columns %in% names(geometry$seeds))) {
        saved_seeds <- as.integer(unlist(
          geometry$seeds[1L, seed_columns, drop = FALSE]
        ))
        expected_seeds <- as.integer(unlist(
          seed_row[1L, seed_columns, drop = FALSE]
        ))

        if (!identical(saved_seeds, expected_seeds)) {
          stop("Saved competing SAR geometry has an incompatible seed record.")
        }
      }
    }

    X_raw <- as.matrix(geometry$X_raw)
    X_cov <- as.matrix(geometry$X_std)
    A_network <- sanitize_network_adjacency(geometry$A_network)

    if (!is.null(geometry$Z_network_modularity)) {
      Z_network <- as.matrix(geometry$Z_network_modularity)
      modularity_q <- ncol(Z_network)
      selected_eigenvalues <- geometry$modularity_selected_eigenvalues
    } else {
      modularity <- modularity_features_single_network(
        A = A_network,
        seed = seed_row$modularity_seed,
        eig_tol = MOD_EIG_TOL,
        max_q = MOD_MAX_Q,
        k_max = MOD_K_MAX
      )
      Z_network <- as.matrix(modularity$Z)
      modularity_q <- modularity$q
      selected_eigenvalues <- modularity$selected_eigenvalues
    }

    sar_selection_summary <- geometry$sar_selection_summary
    sar_selection_diagnostics <- geometry$sar_selection_diagnostics

    if (is.null(sar_selection_summary) ||
        is.null(sar_selection_diagnostics)) {
      reconstructed_sar <- build_sar_label_design(
        X_raw = X_raw,
        sar_seed = seed_row$sar_seed
      )

      if (
        !identical(reconstructed_sar$observed_y, saved_labels$observed_y) ||
        !identical(
          reconstructed_sar$hidden_positive_idx,
          saved_labels$hidden_positive_idx
        )
      ) {
        stop("Saved competing SAR labels cannot be reproduced from X_raw and sar_seed.")
      }

      sar_selection_summary <- reconstructed_sar$selection_summary
      sar_selection_diagnostics <- reconstructed_sar$selection_diagnostics
    }

    return(finalize_feature_bundle(
      source_label = "saved_competing_SAR_geometry",
      X_raw = X_raw,
      X_cov = X_cov,
      A_network = A_network,
      Z_network = Z_network,
      network_summary = geometry$network_summary,
      modularity_q = modularity_q,
      modularity_selected_eigenvalues = selected_eigenvalues,
      seeds = seed_row,
      observed_y = saved_labels$observed_y,
      hidden_positive_idx = saved_labels$hidden_positive_idx,
      observed_positive_idx = saved_labels$observed_positive_idx,
      signal_type = saved_labels$signal_type,
      sar_selection_summary = sar_selection_summary,
      sar_selection_diagnostics = sar_selection_diagnostics,
      covariate_profile = geometry$covariate_profile,
      X_center = geometry$X_center,
      X_scale = geometry$X_scale,
      calibration_label = calibration_label
    ))
  }

  make_generated_feature_bundle(
    seed_row = seed_row,
    source_label = "regenerated_from_competing_SAR_seeds",
    calibration_label = calibration_label
  )
}

## ============================================================================
## 12. FIVE INDEPENDENT GLOBAL PILOTS
## ============================================================================

PRIOR_CALIBRATION_RULES <- list(
  correlation_mode = "relative_to_current_attainable_range",
  corr_pct_range_all = CORR_PCT_RANGE_ALL,
  ell_interval_probs = ELL_INTERVAL_PROBS,
  sigma_interval_probs = SIGMA_INTERVAL_PROBS,
  effective_sd_center = EFFECTIVE_SD_CENTER,
  sigma_range_factor = SIGMA_RANGE_FACTOR,
  ell_bounds_cov = ELL_BOUNDS_COV,
  ell_bounds_network = ELL_BOUNDS_NET,
  n_ell_grid_coarse = N_ELL_GRID_COARSE,
  root_tol = ROOT_TOL,
  half_t_df = HALFT_DF,
  mean_rms_95_mult = MEAN_RMS_95_MULT,
  round_active_priors = ROUND_ACTIVE_PRIORS,
  prior_round_digits = PRIOR_ROUND_DIGITS,
  beta0_prior = beta0_prior_for_sampler
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
    covariate_signal_idx = COVARIATE_SIGNAL_IDX,
    network_signal_idx = NETWORK_SIGNAL_IDX,
    n_hidden_per_pilot = N_HIDDEN,
    n_observed_positive_per_pilot = N_OBSERVED_POS
  ),
  label_mechanism = list(
    type = "covariate-driven SAR",
    exactly_n_hidden_without_replacement = N_HIDDEN,
    lambda = LAMBDA_SAR,
    positive_score_standardization = "within the 45 true positives",
    z_clip = SAR_Z_CLIP,
    log_hiding_weight =
      "-lambda * clipped standardized oracle covariate score",
    oracle_score_available_to_fitted_model = FALSE
  ),
  covariate_dgp = list(
    mu_linear = MU_LINEAR,
    sd_linear = SD_LINEAR,
    block2_x3_offset = BLOCK2_X3_OFFSET,
    block2_x4_center = BLOCK2_X4_CENTER,
    block2_x4_offset = BLOCK2_X4_OFFSET,
    block2_zero_sd_x3 = BLOCK2_ZERO_SD_X3,
    block2_zero_sd_x4 = BLOCK2_ZERO_SD_X4,
    block2_pos_sd_x3 = BLOCK2_POS_SD_X3,
    block2_pos_sd_x4 = BLOCK2_POS_SD_X4,
    disc_radius = DISC_RADIUS,
    ring_radius = RING_RADIUS,
    ring_sd = RING_SD,
    parabola_coef = PARABOLA_COEF,
    delta_parabola = DELTA_PARABOLA,
    sd_parabola_zero = SD_PARABOLA_ZERO,
    sd_parabola_pos = SD_PARABOLA_POS,
    rho_zero = RHO_ZERO,
    rho_pos = RHO_POS
  ),
  network_dgp = list(
    alpha_background = ALPHA_BACKGROUND,
    alpha_signal = ALPHA_FRAUD,
    alpha_legal = ALPHA_LEGAL,
    d_background = D_BACKGROUND,
    d_signal = D_FRAUD,
    d_legal = D_LEGAL,
    triad_closure_prob = TRIAD_CLOSURE_PROB,
    triad_max_new_edges = TRIAD_MAX_NEW_EDGES
  ),
  modularity = list(
    eig_tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q,
    k_max = MOD_K_MAX
  ),
  prior_calibration_rules = PRIOR_CALIBRATION_RULES,
  sampler_controls = list(
    min_blocked_warmup = MIN_BLOCKED_WARMUP,
    max_log_sigma_cov = MAX_LOG_SIGMA_COV,
    max_mean_sd_network = MAX_MEAN_SD_NETWORK,
    max_log_sigma_network = MAX_LOG_SIGMA_NETWORK,
    min_log_sigma = MIN_LOG_SIGMA,
    logl_cov_bounds = LOGL_COV_BOUNDS,
    logl_network_bounds = LOGL_NETWORK_BOUNDS,
    targets = c(f = TARGET_F, theta = TARGET_THETA, joint = TARGET_JOINT)
  )
)

run_or_load_one_pilot <- function(pilot_number) {
  checkpoint_file <- file.path(
    PILOT_CHECKPOINT_DIR,
    sprintf("PoE_pilot_%02d.rds", pilot_number)
  )

  if (
    isTRUE(RESUME_PILOTS) &&
    !isTRUE(FORCE_FRESH_PILOTS) &&
    file.exists(checkpoint_file)
  ) {
    candidate <- tryCatch(
      readRDS(checkpoint_file),
      error = function(e) NULL
    )

    if (
      is.list(candidate) &&
      identical(candidate$settings_signature, PILOT_SETTINGS_SIGNATURE) &&
      identical(as.integer(candidate$pilot), as.integer(pilot_number)) &&
      isTRUE(candidate$complete)
    ) {
      cat(sprintf(
        "[%s] Reusing self-calibrated PoE SAR pilot %02d / %02d\n",
        timestamp_now(),
        pilot_number,
        N_PILOT
      ))
      return(candidate)
    }
  }

  seed_row <- pilot_seed_table[
    pilot_seed_table$pilot == pilot_number,
    ,
    drop = FALSE
  ]

  if (nrow(seed_row) != 1L) {
    stop("Pilot seed lookup failed.")
  }

  bundle <- make_generated_feature_bundle(
    seed_row = seed_row,
    source_label = "independent_global_SAR_pilot",
    calibration_label = sprintf("pilot_%02d", pilot_number)
  )

  pilot_prior_summary <- bundle$prior_calibration$summary
  pilot_prior_summary$stage <- "pilot"
  pilot_prior_summary$index <- pilot_number
  pilot_prior_summary <- pilot_prior_summary[
    ,
    c(
      "stage",
      "index",
      setdiff(names(pilot_prior_summary), c("stage", "index"))
    ),
    drop = FALSE
  ]

  write.csv(
    pilot_prior_summary,
    file.path(
      PRIOR_CALIBRATION_DIR,
      sprintf("PoE_pilot_%02d_prior_calibration.csv", pilot_number)
    ),
    row.names = FALSE
  )

  if (isTRUE(SAVE_PILOT_GEOMETRIES)) {
    saveRDS(
      list(
        pilot = pilot_number,
        seeds = seed_row,
        T = TRUE_T,
        Y = bundle$observed_y,
        hidden_positive_idx = bundle$hidden_positive_idx,
        observed_positive_idx = bundle$observed_positive_idx,
        positive_signal_group = POSITIVE_SIGNAL_GROUP,
        signal_type = bundle$signal_type,
        X_raw = bundle$X_raw,
        X_cov = bundle$X_cov,
        A_network = bundle$A_network,
        Z_network = bundle$Z_network,
        network_summary = bundle$network_summary,
        modularity_q = bundle$modularity_q,
        modularity_selected_eigenvalues =
          bundle$modularity_selected_eigenvalues,
        covariate_profile = bundle$covariate_profile,
        sar_selection_summary = bundle$sar_selection_summary,
        sar_selection_diagnostics = bundle$sar_selection_diagnostics,
        prior_calibration = bundle$prior_calibration
      ),
      file.path(
        PILOT_GEOMETRY_DIR,
        sprintf("PoE_pilot_%02d_geometry.rds", pilot_number)
      )
    )
  }

  mcmc_seed <- as.integer(
    PILOT_MCMC_SEED_BASE + 10000L * pilot_number
  )

  cat("\n\n####################################################################\n")
  cat(sprintf(
    "# SELF-CALIBRATED POE SAR GLOBAL PILOT %d / %d\n",
    pilot_number,
    N_PILOT
  ))
  cat("####################################################################\n")
  cat("Covariate seed       :", seed_row$covariate_seed, "\n")
  cat("Network seed         :", seed_row$network_seed, "\n")
  cat("Modularity seed      :", seed_row$modularity_seed, "\n")
  cat("SAR seed             :", seed_row$sar_seed, "\n")
  cat(
    "SAR hidden groups    : cov=",
    bundle$sar_selection_summary$hidden_covariates_only,
    ", net=",
    bundle$sar_selection_summary$hidden_network_only,
    ", both=",
    bundle$sar_selection_summary$hidden_both,
    "\n",
    sep = ""
  )
  cat("Network mean columns :", ncol(bundle$Z_network), "\n")
  cat("MCMC seed            :", mcmc_seed, "\n")
  cat("Active kernel priors :\n")
  print(bundle$prior_calibration$active_hyperpriors)
  cat("Active tau scales    :\n")
  print(bundle$prior_calibration$tau_scale)

  set.seed(mcmc_seed)

  fit <- run_full_poe_exact_chain(
    Y = bundle$observed_y,
    X_cov = bundle$X_cov,
    A_network = bundle$A_network,
    Z_network = bundle$Z_network,
    network_geometry = bundle$network_geometry,
    hyperpriors = bundle$prior_calibration$active_hyperpriors,
    m_beta = beta0_prior_for_sampler$mean,
    s_beta = beta0_prior_for_sampler$sd,
    tau2_default = bundle$prior_calibration$tau2_default,
    scale_tau = bundle$prior_calibration$tau_scale,
    n_burnin = PILOT_BURNIN,
    n_sample = 0L,
    thin = 1L,
    blocked_frac = PILOT_BLOCKED_FRAC,
    min_blocked_warmup = MIN_BLOCKED_WARMUP,
    init_state = NULL,
    delta_init = 3.5,
    kappa_cov_init = KAPPA_COV_INIT,
    kappa_network_init = KAPPA_NETWORK_INIT,
    proposal_cor_init = PROPOSAL_COR_INIT,
    adapt_block_cor = ADAPT_BLOCK_CORRELATIONS_PILOTS,
    cor_adapt_start = COR_ADAPT_START,
    cor_adapt_interval = COR_ADAPT_INTERVAL,
    cor_adapt_window = COR_ADAPT_WINDOW,
    cor_shrinkage = COR_SHRINKAGE,
    cor_max_abs = COR_MAX_ABS,
    joint_cov_base_init = NULL,
    joint_theta_scale_init = JOINT_THETA_SCALE_INIT_PILOT,
    joint_cor_shrinkage = JOINT_COR_SHRINKAGE,
    joint_cor_max_abs = JOINT_COR_MAX_ABS,
    f_refresh_after_joint = F_REFRESH_AFTER_JOINT,
    nu_tau = NU_TAU,
    sample_tau_groups = SAMPLE_TAU_GROUPS,
    adapt = ADAPT,
    adapt_start = ADAPT_START,
    adapt_interval = ADAPT_INTERVAL,
    target_f = TARGET_F,
    target_theta = TARGET_THETA,
    target_joint = TARGET_JOINT,
    adapt_rate_delta = ADAPT_RATE_DELTA,
    adapt_rate_theta = ADAPT_RATE_THETA,
    adapt_rate_joint_theta = ADAPT_RATE_JOINT_THETA,
    t0_adapt = T0_ADAPT,
    power_adapt = POWER_ADAPT,
    delta_min = DELTA_MIN,
    delta_max = DELTA_MAX,
    kappa_min = KAPPA_MIN,
    kappa_max = KAPPA_MAX,
    joint_theta_scale_min = JOINT_THETA_SCALE_MIN,
    joint_theta_scale_max = JOINT_THETA_SCALE_MAX,
    emp_cov_tail = EMP_COV_TAIL,
    cov_jitter = COV_JITTER,
    max_log_sigma_cov = MAX_LOG_SIGMA_COV,
    max_mean_sd_network = MAX_MEAN_SD_NETWORK,
    max_log_sigma_network = MAX_LOG_SIGMA_NETWORK,
    min_log_sigma = MIN_LOG_SIGMA,
    logl_cov_bounds = LOGL_COV_BOUNDS,
    logl_network_bounds = LOGL_NETWORK_BOUNDS,
    pi_spike = PI_SPIKE,
    spike_log_sig_mean = SPIKE_LOG_SIG_MEAN,
    spike_log_sig_sd = SPIKE_LOG_SIG_SD,
    jitter = JITTER,
    state_tail_length = PILOT_TAIL_LENGTH,
    predictive_label_seed = NULL,
    phase_label = paste0("PoE SAR pilot ", sprintf("%02d", pilot_number)),
    burnin_progress_every = PILOT_PROGRESS_EVERY,
    sample_progress_every = 0L,
    verbose = TRUE
  )

  result <- list(
    settings_signature = PILOT_SETTINGS_SIGNATURE,
    pilot = pilot_number,
    complete = TRUE,
    geometry_source = bundle$source,
    covariate_seed = seed_row$covariate_seed,
    network_seed = seed_row$network_seed,
    modularity_seed = seed_row$modularity_seed,
    sar_seed = seed_row$sar_seed,
    label_design = list(
      observed_y = bundle$observed_y,
      hidden_positive_idx = bundle$hidden_positive_idx,
      observed_positive_idx = bundle$observed_positive_idx,
      signal_type = bundle$signal_type
    ),
    sar_selection_summary = bundle$sar_selection_summary,
    sar_selection_diagnostics = bundle$sar_selection_diagnostics,
    modularity_q = bundle$modularity_q,
    mcmc_seed = mcmc_seed,
    active_priors = list(
      hyperpriors = bundle$prior_calibration$active_hyperpriors,
      hyperpriors_unrounded =
        bundle$prior_calibration$hyperpriors_unrounded,
      tau_scale = bundle$prior_calibration$tau_scale,
      tau_scale_unrounded =
        bundle$prior_calibration$tau_scale_unrounded,
      tau2_default = bundle$prior_calibration$tau2_default,
      theta_prior_mean = bundle$prior_calibration$theta_prior_mean,
      theta_prior_sd = bundle$prior_calibration$theta_prior_sd
    ),
    prior_calibration_summary = pilot_prior_summary,
    tail = fit$tail,
    final_state = fit$state,
    final_proposal = fit$proposal,
    acceptance = fit$acceptance,
    network_component_summary = fit$network_component_summary,
    runtime_sec = fit$time_fit_sec
  )

  saveRDS(result, checkpoint_file)
  rm(bundle, fit)
  invisible(gc())
  result
}

run_all_pilots <- function() {
  pilot_results <- vector("list", N_PILOT)
  for (pilot_number in seq_len(N_PILOT)) {
    pilot_results[[pilot_number]] <- run_or_load_one_pilot(pilot_number)
  }
  pilot_results
}

build_global_calibration <- function(pilot_results) {
  expected_rows <- N_PILOT * PILOT_TAIL_LENGTH
  theta_names <- c(
    "log_ell_cov",
    "log_sigma_cov",
    "log_ell_network",
    "log_sigma_network"
  )
  group_names <- c("cov", "network")
  
  pooled_theta_raw <- do.call(
    rbind,
    lapply(pilot_results, function(x) x$tail$theta_log)
  )
  pooled_theta_standardized <- do.call(
    rbind,
    lapply(pilot_results, function(x) {
      centered <- sweep(
        x$tail$theta_log,
        2L,
        x$active_priors$theta_prior_mean[theta_names],
        "-"
      )
      sweep(
        centered,
        2L,
        x$active_priors$theta_prior_sd[theta_names],
        "/"
      )
    })
  )
  colnames(pooled_theta_standardized) <- theta_names
  
  pooled_beta_fixed <- do.call(
    rbind,
    lapply(pilot_results, function(x) x$tail$beta_fixed)
  )
  pooled_tau2_ratio <- do.call(
    rbind,
    lapply(pilot_results, function(x) {
      sweep(
        x$tail$tau2[, group_names, drop = FALSE],
        2L,
        x$active_priors$tau_scale[group_names]^2,
        "/"
      )
    })
  )
  pooled_aux_ratio <- do.call(
    rbind,
    lapply(pilot_results, function(x) {
      sweep(
        x$tail$aux_tau[, group_names, drop = FALSE],
        2L,
        x$active_priors$tau_scale[group_names]^2,
        "*"
      )
    })
  )
  
  if (nrow(pooled_theta_raw) != expected_rows ||
      nrow(pooled_theta_standardized) != expected_rows ||
      nrow(pooled_beta_fixed) != expected_rows ||
      nrow(pooled_tau2_ratio) != expected_rows ||
      nrow(pooled_aux_ratio) != expected_rows) {
    stop("The pooled pilot tails have an unexpected number of rows.")
  }
  
  pooled_joint_correlation <- estimate_correlation_matrix(
    pooled_theta_standardized,
    fallback = diag(4L),
    shrinkage = JOINT_COR_SHRINKAGE,
    max_abs = JOINT_COR_MAX_ABS,
    min_n = 200L
  )
  dimnames(pooled_joint_correlation) <- list(theta_names, theta_names)
  
  pooled_block_correlations <- list(
    cov = estimate_cor_2x2(
      pooled_theta_standardized[, 1:2, drop = FALSE],
      fallback = diag(2),
      shrinkage = COR_SHRINKAGE,
      max_abs = COR_MAX_ABS,
      min_n = 200L
    ),
    network = estimate_cor_2x2(
      pooled_theta_standardized[, 3:4, drop = FALSE],
      fallback = diag(2),
      shrinkage = COR_SHRINKAGE,
      max_abs = COR_MAX_ABS,
      min_n = 200L
    )
  )
  
  standardized_final_joint_covariances <- lapply(
    pilot_results,
    function(x) {
      prior_sd <- x$active_priors$theta_prior_sd[theta_names]
      D_inverse <- diag(1 / prior_sd, 4L)
      sanitize_covariance_matrix(
        D_inverse %*%
          x$final_proposal$Sigma_theta_final %*%
          D_inverse
      )
    }
  )
  global_joint_covariance_standardized <- median_covariance_matrix(
    standardized_final_joint_covariances,
    dimension = 4L
  )
  dimnames(global_joint_covariance_standardized) <- list(
    theta_names,
    theta_names
  )
  
  global_delta <- median(vapply(
    pilot_results,
    function(x) x$final_proposal$delta_final,
    numeric(1)
  ))
  global_kappa_cov <- median(vapply(
    pilot_results,
    function(x) x$final_proposal$kappa_cov_final,
    numeric(1)
  ))
  global_kappa_network <- median(vapply(
    pilot_results,
    function(x) x$final_proposal$kappa_network_final,
    numeric(1)
  ))
  
  global_initial_values <- list(
    theta_prior_z_offset = setNames(
      apply(pooled_theta_standardized, 2L, median),
      theta_names
    ),
    beta0 = median(pooled_beta_fixed[, "beta0"]),
    cov_beta = setNames(
      apply(pooled_beta_fixed[, 2:11, drop = FALSE], 2L, median),
      colnames(pooled_beta_fixed)[2:11]
    ),
    tau2_scale_ratio = setNames(
      apply(pooled_tau2_ratio, 2L, median),
      group_names
    ),
    aux_tau_scale_ratio = setNames(
      apply(pooled_aux_ratio, 2L, median),
      group_names
    ),
    delta_init = global_delta,
    kappa_cov_init = global_kappa_cov,
    kappa_network_init = global_kappa_network,
    proposal_cor_init = pooled_block_correlations,
    joint_covariance_standardized =
      global_joint_covariance_standardized,
    joint_theta_scale_init = JOINT_THETA_SCALE_INIT_PRODUCTION,
    pooled_joint_correlation = pooled_joint_correlation,
    pooled_tail_draws = nrow(pooled_theta_standardized)
  )
  
  parameter_vectors <- list(
    z_log_ell_cov = pooled_theta_standardized[, "log_ell_cov"],
    z_log_sigma_cov = pooled_theta_standardized[, "log_sigma_cov"],
    z_log_ell_network = pooled_theta_standardized[, "log_ell_network"],
    z_log_sigma_network = pooled_theta_standardized[, "log_sigma_network"],
    raw_log_ell_cov = pooled_theta_raw[, "log_ell_cov"],
    raw_log_sigma_cov = pooled_theta_raw[, "log_sigma_cov"],
    raw_log_ell_network = pooled_theta_raw[, "log_ell_network"],
    raw_log_sigma_network = pooled_theta_raw[, "log_sigma_network"],
    beta0 = pooled_beta_fixed[, "beta0"],
    tau2_cov_scale_ratio = pooled_tau2_ratio[, "cov"],
    tau2_network_scale_ratio = pooled_tau2_ratio[, "network"],
    aux_cov_scale_ratio = pooled_aux_ratio[, "cov"],
    aux_network_scale_ratio = pooled_aux_ratio[, "network"]
  )
  for (column_name in colnames(pooled_beta_fixed)[2:11]) {
    parameter_vectors[[column_name]] <- pooled_beta_fixed[, column_name]
  }
  
  distribution_rows <- lapply(
    names(parameter_vectors),
    function(parameter_name) {
      values <- parameter_vectors[[parameter_name]]
      statistics <- summary_stats(values)
      data.frame(
        parameter = parameter_name,
        mean = unname(statistics["mean"]),
        sd = unname(statistics["sd"]),
        q025 = unname(statistics["q025"]),
        median = unname(statistics["median"]),
        q975 = unname(statistics["q975"]),
        pooled_tail_draws = length(values),
        stringsAsFactors = FALSE
      )
    }
  )
  
  per_pilot_rows <- lapply(seq_along(pilot_results), function(pp) {
    x <- pilot_results[[pp]]
    hp <- x$active_priors$hyperpriors
    data.frame(
      pilot = pp,
      covariate_seed = x$covariate_seed,
      network_seed = x$network_seed,
      modularity_seed = x$modularity_seed,
      modularity_q = x$modularity_q,
      mcmc_seed = x$mcmc_seed,
      prior_log_ell_cov_mean = hp$expert0_cov$log_ell_mean,
      prior_log_ell_cov_sd = hp$expert0_cov$log_ell_sd,
      prior_log_sigma_cov_mean = hp$expert0_cov$log_sig_mean,
      prior_log_sigma_cov_sd = hp$expert0_cov$log_sig_sd,
      prior_log_ell_network_mean = hp$expert1_network$log_ell_mean,
      prior_log_ell_network_sd = hp$expert1_network$log_ell_sd,
      prior_log_sigma_network_mean = hp$expert1_network$log_sig_mean,
      prior_log_sigma_network_sd = hp$expert1_network$log_sig_sd,
      tau_scale_cov = x$active_priors$tau_scale["cov"],
      tau_scale_network = x$active_priors$tau_scale["network"],
      log_ell_cov_median = median(
        x$tail$theta_log[, "log_ell_cov"]
      ),
      log_sigma_cov_median = median(
        x$tail$theta_log[, "log_sigma_cov"]
      ),
      log_ell_network_median = median(
        x$tail$theta_log[, "log_ell_network"]
      ),
      log_sigma_network_median = median(
        x$tail$theta_log[, "log_sigma_network"]
      ),
      beta0_median = median(x$tail$beta_fixed[, "beta0"]),
      tau2_cov_ratio_median = median(
        x$tail$tau2[, "cov"] /
          x$active_priors$tau_scale["cov"]^2
      ),
      tau2_network_ratio_median = median(
        x$tail$tau2[, "network"] /
          x$active_priors$tau_scale["network"]^2
      ),
      delta_final = x$final_proposal$delta_final,
      kappa_cov_final = x$final_proposal$kappa_cov_final,
      kappa_network_final = x$final_proposal$kappa_network_final,
      rho_cov_final =
        x$final_proposal$proposal_cor_final$cov[1L, 2L],
      rho_network_final =
        x$final_proposal$proposal_cor_final$network[1L, 2L],
      joint_theta_scale_final =
        x$final_proposal$joint_theta_scale_final,
      blocked_acceptance_f = x$acceptance$blocked_f,
      blocked_acceptance_cov = x$acceptance$blocked_cov,
      blocked_acceptance_network = x$acceptance$blocked_network,
      joint_tune_acceptance = x$acceptance$joint_tune,
      runtime_sec = x$runtime_sec,
      stringsAsFactors = FALSE
    )
  })
  
  pilot_prior_summary <- do.call(
    rbind,
    lapply(pilot_results, `[[`, "prior_calibration_summary")
  )
  row.names(pilot_prior_summary) <- NULL
  
  list(
    settings_signature = PILOT_SETTINGS_SIGNATURE,
    created_at = timestamp_now(),
    global_initial_values = global_initial_values,
    distribution_summary = do.call(rbind, distribution_rows),
    per_pilot_summary = do.call(rbind, per_pilot_rows),
    pilot_prior_calibration_summary = pilot_prior_summary,
    pilot_seed_table = pilot_seed_table
  )
}

load_or_build_global_calibration <- function() {
  if (isTRUE(RESUME_PILOTS) && !isTRUE(FORCE_FRESH_PILOTS) &&
      file.exists(GLOBAL_CALIBRATION_FILE)) {
    candidate <- tryCatch(
      readRDS(GLOBAL_CALIBRATION_FILE),
      error = function(e) NULL
    )
    if (is.list(candidate) &&
        identical(candidate$settings_signature, PILOT_SETTINGS_SIGNATURE)) {
      cat("\nLoaded the completed five-pilot self-calibrated PoE run.\n")
      return(candidate)
    }
  }
  
  pilot_results <- run_all_pilots()
  calibration <- build_global_calibration(pilot_results)
  saveRDS(calibration, GLOBAL_CALIBRATION_FILE)
  
  global_pilot_calibration <- calibration
  save(global_pilot_calibration, file = GLOBAL_CALIBRATION_RDATA)
  
  write.csv(
    calibration$distribution_summary,
    file.path(OUT_DIR, "PoE_global_pilot_distribution_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    calibration$per_pilot_summary,
    file.path(OUT_DIR, "PoE_global_pilot_per_pilot_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    calibration$pilot_prior_calibration_summary,
    file.path(OUT_DIR, "PoE_pilot_prior_calibration_summary.csv"),
    row.names = FALSE
  )
  
  rm(pilot_results)
  invisible(gc())
  calibration
}

global_pilot_calibration <- load_or_build_global_calibration()
GLOBAL_CALIBRATION_HASH <- file_md5_or_na(GLOBAL_CALIBRATION_FILE)

write.csv(
  global_pilot_calibration$distribution_summary,
  file.path(OUT_DIR, "PoE_global_pilot_distribution_summary.csv"),
  row.names = FALSE
)
write.csv(
  global_pilot_calibration$per_pilot_summary,
  file.path(OUT_DIR, "PoE_global_pilot_per_pilot_summary.csv"),
  row.names = FALSE
)
write.csv(
  global_pilot_calibration$pilot_prior_calibration_summary,
  file.path(OUT_DIR, "PoE_pilot_prior_calibration_summary.csv"),
  row.names = FALSE
)

cat("\n================ GLOBAL PRIOR-RELATIVE INITIAL VALUES ================\n")
print(global_pilot_calibration$global_initial_values)

## ============================================================================
## 13. GEOMETRY-SPECIFIC PRODUCTION INITIALIZATION AND EVALUATION
## ============================================================================

make_global_initialization <- function(bundle, calibration) {
  q_network <- ncol(bundle$Z_network)
  if (q_network < 1L) {
    stop("The current production network has no modularity features.")
  }
  
  prior <- bundle$prior_calibration
  theta_names <- names(prior$theta_prior_mean)
  group_names <- c("cov", "network")
  
  theta_initial <- prior$theta_prior_mean[theta_names] +
    prior$theta_prior_sd[theta_names] *
    calibration$theta_prior_z_offset[theta_names]
  
  tau2_initial <- prior$tau_scale[group_names]^2 *
    calibration$tau2_scale_ratio[group_names]
  aux_initial <- calibration$aux_tau_scale_ratio[group_names] /
    prior$tau_scale[group_names]^2
  
  beta <- c(
    calibration$beta0,
    as.numeric(calibration$cov_beta),
    rep(0, q_network)
  )
  
  if (length(beta) != 1L + 10L + q_network ||
      any(!is.finite(beta)) ||
      any(!is.finite(theta_initial)) ||
      any(!is.finite(tau2_initial)) || any(tau2_initial <= 0) ||
      any(!is.finite(aux_initial)) || any(aux_initial <= 0)) {
    stop("The geometry-specific global initialization is invalid.")
  }
  
  D_prior <- diag(prior$theta_prior_sd[theta_names], 4L)
  joint_covariance_current <- sanitize_covariance_matrix(
    D_prior %*%
      calibration$joint_covariance_standardized %*%
      D_prior
  )
  dimnames(joint_covariance_current) <- list(theta_names, theta_names)
  
  list(
    state = list(
      theta_log = setNames(as.numeric(theta_initial), theta_names),
      beta = beta,
      tau2 = setNames(as.numeric(tau2_initial), group_names),
      aux_tau = setNames(as.numeric(aux_initial), group_names),
      ## No current-data latent f is transferred between data sets.
      f = NULL
    ),
    joint_cov_base_init = joint_covariance_current
  )
}

evaluate_full_poe <- function(
    fit,
    replicate_number,
    runtime_sec,
    observed_y,
    hidden_positive_idx,
    signal_type
) {
  labels <- validate_sar_label_bundle(
    observed_y = observed_y,
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = which(observed_y == 1L),
    signal_type = signal_type
  )

  observed_y <- labels$observed_y
  hidden_positive_idx <- labels$hidden_positive_idx
  signal_type <- labels$signal_type

  probabilities <- validate_probability_vector(
    fit$posterior_mean_scores,
    N_NODES,
    METHOD_NAME
  )

  predictive_T_draws <- as.matrix(fit$predictive_T_draws)

  if (!all(dim(predictive_T_draws) == c(fit$n_saved, N_NODES))) {
    stop("The predictive label draws have incompatible dimensions.")
  }

  if (any(!predictive_T_draws %in% c(0L, 1L))) {
    stop("The predictive label draws must be binary.")
  }

  unlabeled_idx <- which(observed_y == 0L)
  hidden_truth_unlabeled <- as.integer(
    TRUE_T[unlabeled_idx] == 1L
  )

  if (
    length(unlabeled_idx) != N_UNLABELED ||
    sum(hidden_truth_unlabeled) != N_HIDDEN ||
    !identical(
      sort(unlabeled_idx[hidden_truth_unlabeled == 1L]),
      hidden_positive_idx
    )
  ) {
    stop("The replication-specific SAR unlabeled evaluation set is inconsistent.")
  }

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

  overall_80 <- predictive_set_metrics(
    predictive_T_draws = predictive_T_draws,
    truth = TRUE_T,
    evaluation_idx = seq_len(N_NODES),
    level = 0.80
  )

  hidden_80 <- predictive_set_metrics(
    predictive_T_draws = predictive_T_draws,
    truth = TRUE_T,
    evaluation_idx = hidden_positive_idx,
    level = 0.80
  )

  overall_95 <- predictive_set_metrics(
    predictive_T_draws = predictive_T_draws,
    truth = TRUE_T,
    evaluation_idx = seq_len(N_NODES),
    level = 0.95
  )

  hidden_95 <- predictive_set_metrics(
    predictive_T_draws = predictive_T_draws,
    truth = TRUE_T,
    evaluation_idx = hidden_positive_idx,
    level = 0.95
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

    Hidden_LogLoss = overall_logloss(
      probabilities[hidden_positive_idx],
      TRUE_T[hidden_positive_idx]
    ),

    Overall_Brier = overall_brier(
      probabilities,
      TRUE_T
    ),

    Hidden_Brier = overall_brier(
      probabilities[hidden_positive_idx],
      TRUE_T[hidden_positive_idx]
    ),

    Overall_Coverage_80 = overall_80$coverage,
    Hidden_Coverage_80 = hidden_80$coverage,
    Overall_Interval_Score_80 = overall_80$mean_interval_score,
    Hidden_Interval_Score_80 = hidden_80$mean_interval_score,
    Overall_Interval_Width_80 = overall_80$mean_width,
    Hidden_Interval_Width_80 = hidden_80$mean_width,

    Overall_Coverage_95 = overall_95$coverage,
    Hidden_Coverage_95 = hidden_95$coverage,
    Overall_Interval_Score_95 = overall_95$mean_interval_score,
    Hidden_Interval_Score_95 = hidden_95$mean_interval_score,
    Overall_Interval_Width_95 = overall_95$mean_width,
    Hidden_Interval_Width_95 = hidden_95$mean_width,

    anchor_max_mean = mean(fit$max_probability_draws),
    anchor_max_median = median(fit$max_probability_draws),
    anchor_max_q025 = unname(
      quantile(fit$max_probability_draws, 0.025)
    ),
    anchor_max_q975 = unname(
      quantile(fit$max_probability_draws, 0.975)
    ),
    anchor_Pr_over_090 = mean(fit$max_probability_draws > 0.90),
    anchor_Pr_over_095 = mean(fit$max_probability_draws > 0.95),
    anchor_Pr_over_099 = mean(fit$max_probability_draws > 0.99),
    anchor_q99_median = median(fit$q99_probability_draws),
    anchor_top5_median = median(fit$top5_mean_probability_draws),

    beta0_posterior_mean = mean(fit$beta0_draws),
    tau_cov_posterior_mean = mean(fit$tau_draws[, "tau_cov"]),
    tau_network_posterior_mean = mean(
      fit$tau_draws[, "tau_network"]
    ),
    ell_cov_posterior_mean = mean(fit$theta_draws[, "ell_cov"]),
    sigma_cov_posterior_mean = mean(fit$theta_draws[, "sigma_cov"]),
    ell_network_posterior_mean = mean(
      fit$theta_draws[, "ell_network"]
    ),
    sigma_network_posterior_mean = mean(
      fit$theta_draws[, "sigma_network"]
    ),
    network_mean_feature_count = fit$network_feature_count,

    acc_blocked_f = fit$acceptance$blocked_f,
    acc_blocked_cov = fit$acceptance$blocked_cov,
    acc_blocked_network = fit$acceptance$blocked_network,
    acc_joint_tune = fit$acceptance$joint_tune,
    acc_joint_sampling = fit$acceptance$joint_sampling,
    acc_f_refresh_sampling = fit$acceptance$f_refresh_sampling,
    delta_final = fit$proposal$delta_final,
    kappa_cov_final = fit$proposal$kappa_cov_final,
    kappa_network_final = fit$proposal$kappa_network_final,
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
    positive_signal_group = POSITIVE_SIGNAL_GROUP,
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

    PI80_lower = overall_80$lower,
    PI80_upper = overall_80$upper,
    PI80_width = overall_80$upper - overall_80$lower,
    covered_80 = overall_80$covered,
    interval_score_80 = overall_80$interval_score,

    PI95_lower = overall_95$lower,
    PI95_upper = overall_95$upper,
    PI95_width = overall_95$upper - overall_95$lower,
    covered_95 = overall_95$covered,
    interval_score_95 = overall_95$interval_score,

    stringsAsFactors = FALSE
  )

  node_scores$unlabeled_rank[unlabeled_idx] <- unlabeled_rank

  list(
    metrics = metrics,
    node_scores = node_scores
  )
}

## ============================================================================
## 14. PRODUCTION SETTINGS SIGNATURE
## ============================================================================

PRODUCTION_SETTINGS_SIGNATURE <- list(
  script_version = SCRIPT_VERSION,
  master_seed = MASTER_SEED,
  n_rep = N_REP,
  production_burnin = PRODUCTION_BURNIN,
  production_sampling = PRODUCTION_SAMPLING,
  production_thin = PRODUCTION_THIN,
  production_blocked_frac = PRODUCTION_BLOCKED_FRAC,
  global_calibration_hash = GLOBAL_CALIBRATION_HASH,
  use_saved_competing_geometries = USE_SAVED_COMPETING_GEOMETRIES,
  data_design = list(
    true_t = TRUE_T,
    true_zero_idx = TRUE_ZERO_IDX,
    true_positive_idx = TRUE_POS_IDX,
    covariates_only_idx = COV_ONLY_IDX,
    network_only_idx = NET_ONLY_IDX,
    both_idx = BOTH_IDX,
    covariate_signal_idx = COVARIATE_SIGNAL_IDX,
    network_signal_idx = NETWORK_SIGNAL_IDX,
    n_hidden_per_replication = N_HIDDEN,
    n_observed_positive_per_replication = N_OBSERVED_POS
  ),
  label_mechanism = list(
    type = "covariate-driven SAR",
    exactly_n_hidden_without_replacement = N_HIDDEN,
    lambda = LAMBDA_SAR,
    positive_score_standardization = "within the 45 true positives",
    z_clip = SAR_Z_CLIP,
    log_hiding_weight =
      "-lambda * clipped standardized oracle covariate score",
    oracle_score_available_to_fitted_model = FALSE
  ),
  covariate_dgp = PILOT_SETTINGS_SIGNATURE$covariate_dgp,
  network_dgp = PILOT_SETTINGS_SIGNATURE$network_dgp,
  modularity = PILOT_SETTINGS_SIGNATURE$modularity,
  prior_calibration_rules = PRIOR_CALIBRATION_RULES,
  sampler_controls = list(
    adapt_block_cor = ADAPT_BLOCK_CORRELATIONS_PRODUCTION,
    targets = c(f = TARGET_F, theta = TARGET_THETA, joint = TARGET_JOINT),
    max_log_sigma_cov = MAX_LOG_SIGMA_COV,
    max_mean_sd_network = MAX_MEAN_SD_NETWORK,
    max_log_sigma_network = MAX_LOG_SIGMA_NETWORK,
    min_log_sigma = MIN_LOG_SIGMA,
    logl_cov_bounds = LOGL_COV_BOUNDS,
    logl_network_bounds = LOGL_NETWORK_BOUNDS
  ),
  evaluation = list(
    version = EVALUATION_VERSION,
    hidden_top_k = HIDDEN_TOP_K,
    overall_recall_k = OVERALL_RECALL_K,
    predictive_interval_levels = PREDICTIVE_INTERVAL_LEVELS,
    predictive_label_distribution = paste0(
      "For every node and retained posterior state, draw ",
      "T_pred ~ Bernoulli(plogis(f)) under a separate predictive seed."
    ),
    interval_quantiles = "central discrete empirical quantiles with type=1",
    ground_truth = paste0(
      "fixed TRUE_T; every SAR-hidden positive and every observed positive ",
      "has T=1"
    ),
    H_intervals = FALSE,
    poisson_binomial = FALSE
  )
)

## ============================================================================
## 15. ONE PRODUCTION REPLICATION
## ============================================================================

run_one_production_replication <- function(replicate_number) {
  checkpoint_file <- file.path(
    PRODUCTION_CHECKPOINT_DIR,
    sprintf("PoE_replicate_%03d.rds", replicate_number)
  )

  if (
    isTRUE(RESUME_PRODUCTION) &&
    !isTRUE(FORCE_FRESH_PRODUCTION) &&
    file.exists(checkpoint_file)
  ) {
    candidate <- tryCatch(
      readRDS(checkpoint_file),
      error = function(e) NULL
    )

    if (
      is.list(candidate) &&
      identical(candidate$settings_signature, PRODUCTION_SETTINGS_SIGNATURE) &&
      identical(as.integer(candidate$replicate), as.integer(replicate_number)) &&
      isTRUE(candidate$complete)
    ) {
      cat(sprintf(
        "[%s] Reusing self-calibrated PoE SAR replication %03d / %03d\n",
        timestamp_now(),
        replicate_number,
        N_REP
      ))
      return(candidate)
    }
  }

  replication_start <- Sys.time()

  bundle <- load_or_generate_production_bundle(replicate_number)

  sar_selection_summary_now <- bundle$sar_selection_summary
  sar_selection_summary_now$replicate <- replicate_number
  sar_selection_summary_now <- sar_selection_summary_now[
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

  sar_selection_diagnostics_now <- bundle$sar_selection_diagnostics
  sar_selection_diagnostics_now$replicate <- replicate_number
  sar_selection_diagnostics_now <- sar_selection_diagnostics_now[
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

  global_values <- global_pilot_calibration$global_initial_values
  initialization <- make_global_initialization(bundle, global_values)
  init_state <- initialization$state

  mcmc_seed <- as.integer(
    PRODUCTION_MCMC_SEED_BASE + 10000L * replicate_number
  )

  predictive_seed <- as.integer(
    PRODUCTION_PREDICTIVE_SEED_BASE + 10000L * replicate_number
  )

  cat("\n\n####################################################################\n")
  cat(sprintf(
    "# SELF-CALIBRATED FULL POE SAR REPLICATION %d / %d\n",
    replicate_number,
    N_REP
  ))
  cat("####################################################################\n")
  cat("Geometry source       :", bundle$source, "\n")
  cat("Covariate seed        :", bundle$seeds$covariate_seed, "\n")
  cat("Network seed          :", bundle$seeds$network_seed, "\n")
  cat("Modularity seed       :", bundle$seeds$modularity_seed, "\n")
  cat("SAR seed              :", bundle$seeds$sar_seed, "\n")
  cat(
    "SAR hidden groups     : cov=",
    sar_selection_summary_now$hidden_covariates_only,
    ", net=",
    sar_selection_summary_now$hidden_network_only,
    ", both=",
    sar_selection_summary_now$hidden_both,
    "\n",
    sep = ""
  )
  cat(
    "SAR mean scores       : hidden=",
    signif(sar_selection_summary_now$mean_score_hidden, 6),
    ", observed=",
    signif(sar_selection_summary_now$mean_score_observed, 6),
    "\n",
    sep = ""
  )
  cat("Network mean features :", ncol(bundle$Z_network), "\n")
  cat(
    "Prior calibration time:",
    format_duration(bundle$prior_calibration$runtime_sec),
    "\n"
  )
  cat("MCMC seed             :", mcmc_seed, "\n")
  cat("Predictive-label seed :", predictive_seed, "\n")
  cat("Active kernel priors  :\n")
  print(bundle$prior_calibration$active_hyperpriors)
  cat("Active tau scales     :\n")
  print(bundle$prior_calibration$tau_scale)

  fit_error <- NA_character_

  fit <- tryCatch({
    set.seed(mcmc_seed)

    run_full_poe_exact_chain(
      Y = bundle$observed_y,
      X_cov = bundle$X_cov,
      A_network = bundle$A_network,
      Z_network = bundle$Z_network,
      network_geometry = bundle$network_geometry,
      hyperpriors = bundle$prior_calibration$active_hyperpriors,
      m_beta = beta0_prior_for_sampler$mean,
      s_beta = beta0_prior_for_sampler$sd,
      tau2_default = bundle$prior_calibration$tau2_default,
      scale_tau = bundle$prior_calibration$tau_scale,
      n_burnin = PRODUCTION_BURNIN,
      n_sample = PRODUCTION_SAMPLING,
      thin = PRODUCTION_THIN,
      blocked_frac = PRODUCTION_BLOCKED_FRAC,
      min_blocked_warmup = MIN_BLOCKED_WARMUP,
      init_state = init_state,
      delta_init = global_values$delta_init,
      kappa_cov_init = global_values$kappa_cov_init,
      kappa_network_init = global_values$kappa_network_init,
      proposal_cor_init = global_values$proposal_cor_init,
      adapt_block_cor = ADAPT_BLOCK_CORRELATIONS_PRODUCTION,
      cor_adapt_start = COR_ADAPT_START,
      cor_adapt_interval = COR_ADAPT_INTERVAL,
      cor_adapt_window = COR_ADAPT_WINDOW,
      cor_shrinkage = COR_SHRINKAGE,
      cor_max_abs = COR_MAX_ABS,
      joint_cov_base_init = initialization$joint_cov_base_init,
      joint_theta_scale_init = global_values$joint_theta_scale_init,
      joint_cor_shrinkage = JOINT_COR_SHRINKAGE,
      joint_cor_max_abs = JOINT_COR_MAX_ABS,
      f_refresh_after_joint = F_REFRESH_AFTER_JOINT,
      nu_tau = NU_TAU,
      sample_tau_groups = SAMPLE_TAU_GROUPS,
      adapt = ADAPT,
      adapt_start = ADAPT_START,
      adapt_interval = ADAPT_INTERVAL,
      target_f = TARGET_F,
      target_theta = TARGET_THETA,
      target_joint = TARGET_JOINT,
      adapt_rate_delta = ADAPT_RATE_DELTA,
      adapt_rate_theta = ADAPT_RATE_THETA,
      adapt_rate_joint_theta = ADAPT_RATE_JOINT_THETA,
      t0_adapt = T0_ADAPT,
      power_adapt = POWER_ADAPT,
      delta_min = DELTA_MIN,
      delta_max = DELTA_MAX,
      kappa_min = KAPPA_MIN,
      kappa_max = KAPPA_MAX,
      joint_theta_scale_min = JOINT_THETA_SCALE_MIN,
      joint_theta_scale_max = JOINT_THETA_SCALE_MAX,
      emp_cov_tail = EMP_COV_TAIL,
      cov_jitter = COV_JITTER,
      max_log_sigma_cov = MAX_LOG_SIGMA_COV,
      max_mean_sd_network = MAX_MEAN_SD_NETWORK,
      max_log_sigma_network = MAX_LOG_SIGMA_NETWORK,
      min_log_sigma = MIN_LOG_SIGMA,
      logl_cov_bounds = LOGL_COV_BOUNDS,
      logl_network_bounds = LOGL_NETWORK_BOUNDS,
      pi_spike = PI_SPIKE,
      spike_log_sig_mean = SPIKE_LOG_SIG_MEAN,
      spike_log_sig_sd = SPIKE_LOG_SIG_SD,
      jitter = JITTER,
      state_tail_length = 0L,
      predictive_label_seed = predictive_seed,
      phase_label = paste0(
        "PoE SAR replication ",
        sprintf("%03d", replicate_number)
      ),
      burnin_progress_every = PRODUCTION_BURNIN_PROGRESS_EVERY,
      sample_progress_every = PRODUCTION_SAMPLE_PROGRESS_EVERY,
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
      settings_signature = PRODUCTION_SETTINGS_SIGNATURE,
      replicate = replicate_number,
      complete = FALSE,
      failure = failure,
      sar_selection_summary = sar_selection_summary_now,
      sar_selection_diagnostics = sar_selection_diagnostics_now,
      label_design = list(
        observed_y = bundle$observed_y,
        hidden_positive_idx = bundle$hidden_positive_idx,
        observed_positive_idx = bundle$observed_positive_idx,
        signal_type = bundle$signal_type
      ),
      prior_calibration_summary = bundle$prior_calibration$summary,
      time_start = replication_start,
      time_end = Sys.time(),
      runtime_sec = runtime_sec
    )

    saveRDS(result, checkpoint_file)

    if (isTRUE(FAIL_IF_PRODUCTION_FAILS)) {
      stop(
        "PoE SAR failed in replication ",
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

  evaluation <- evaluate_full_poe(
    fit = fit,
    replicate_number = replicate_number,
    runtime_sec = runtime_sec,
    observed_y = bundle$observed_y,
    hidden_positive_idx = bundle$hidden_positive_idx,
    signal_type = bundle$signal_type
  )

  prior_summary <- bundle$prior_calibration$summary
  prior_summary$scenario <- SCENARIO_NAME
  prior_summary$replicate <- replicate_number
  prior_summary$geometry_source <- bundle$source
  prior_summary <- prior_summary[
    ,
    c(
      "scenario",
      "replicate",
      "geometry_source",
      setdiff(
        names(prior_summary),
        c("scenario", "replicate", "geometry_source")
      )
    ),
    drop = FALSE
  ]

  write.csv(
    prior_summary,
    file.path(
      PRIOR_CALIBRATION_DIR,
      sprintf("PoE_replicate_%03d_prior_calibration.csv", replicate_number)
    ),
    row.names = FALSE
  )

  network_summary <- bundle$network_summary

  if (is.null(network_summary)) {
    network_summary <- data.frame(stringsAsFactors = FALSE)
  } else {
    network_summary <- as.data.frame(
      network_summary,
      stringsAsFactors = FALSE
    )
  }

  geometry_base <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    geometry_source = bundle$source,
    covariate_seed = bundle$seeds$covariate_seed,
    network_seed = bundle$seeds$network_seed,
    modularity_seed = bundle$seeds$modularity_seed,
    sar_seed = bundle$seeds$sar_seed,
    lambda_SAR = LAMBDA_SAR,
    n = N_NODES,
    n_covariates = ncol(bundle$X_cov),
    n_network_features = ncol(bundle$Z_network),
    n_true_zero = N_TRUE_ZERO,
    n_true_positive = N_TRUE_POS,
    n_observed_positive = length(bundle$observed_positive_idx),
    n_hidden_positive = length(bundle$hidden_positive_idx),
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
    n_covariate_signal = length(COVARIATE_SIGNAL_IDX),
    n_network_signal = length(NETWORK_SIGNAL_IDX),
    n_both_signal = length(intersect(
      COVARIATE_SIGNAL_IDX,
      NETWORK_SIGNAL_IDX
    )),
    prior_calibration_runtime_sec = bundle$prior_calibration$runtime_sec,
    tau_scale_cov = bundle$prior_calibration$tau_scale["cov"],
    tau_scale_network = bundle$prior_calibration$tau_scale["network"],
    stringsAsFactors = FALSE
  )

  geometry_summary <- if (ncol(network_summary) > 0L) {
    cbind(geometry_base, network_summary)
  } else {
    geometry_base
  }

  modularity_summary <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    q_network = bundle$modularity_q,
    Z_network_rms_factor = sqrt(
      sum(bundle$Z_network^2) / nrow(bundle$Z_network)
    ),
    tau_scale_network = bundle$prior_calibration$tau_scale["network"],
    selected_eigenvalues = if (
      length(bundle$modularity_selected_eigenvalues) > 0L
    ) {
      paste(
        signif(bundle$modularity_selected_eigenvalues, 8),
        collapse = ";"
      )
    } else {
      ""
    },
    stringsAsFactors = FALSE
  )

  result <- list(
    settings_signature = PRODUCTION_SETTINGS_SIGNATURE,
    replicate = replicate_number,
    complete = TRUE,
    metrics = evaluation$metrics,
    node_scores = evaluation$node_scores,
    sar_selection_summary = sar_selection_summary_now,
    sar_selection_diagnostics = sar_selection_diagnostics_now,
    label_design = list(
      observed_y = bundle$observed_y,
      hidden_positive_idx = bundle$hidden_positive_idx,
      observed_positive_idx = bundle$observed_positive_idx,
      signal_type = bundle$signal_type
    ),
    prior_calibration_summary = prior_summary,
    active_priors = list(
      hyperpriors = bundle$prior_calibration$active_hyperpriors,
      hyperpriors_unrounded =
        bundle$prior_calibration$hyperpriors_unrounded,
      tau_scale = bundle$prior_calibration$tau_scale,
      tau_scale_unrounded =
        bundle$prior_calibration$tau_scale_unrounded,
      tau2_default = bundle$prior_calibration$tau2_default,
      theta_prior_mean = bundle$prior_calibration$theta_prior_mean,
      theta_prior_sd = bundle$prior_calibration$theta_prior_sd
    ),
    geometry_summary = geometry_summary,
    modularity_summary = modularity_summary,
    final_state = fit$state,
    final_proposal = fit$proposal,
    acceptance = fit$acceptance,
    network_component_summary = fit$network_component_summary,
    n_saved = fit$n_saved,
    mcmc_seed = mcmc_seed,
    predictive_label_seed = predictive_seed,
    time_start = replication_start,
    time_end = Sys.time(),
    runtime_sec = runtime_sec,
    failure = NULL
  )

  saveRDS(result, checkpoint_file)

  rm(bundle, initialization, init_state, fit, evaluation)
  invisible(gc())
  result
}

## ============================================================================
## 16. RUN ALL 100 PRODUCTION REPLICATIONS
## ============================================================================

production_results <- vector("list", N_REP)
production_start <- Sys.time()
n_loaded_complete <- 0L
n_computed_now <- 0L

for (replicate_number in seq_len(N_REP)) {
  checkpoint_file <- file.path(
    PRODUCTION_CHECKPOINT_DIR,
    sprintf("PoE_replicate_%03d.rds", replicate_number)
  )
  
  was_complete_before <- FALSE
  if (isTRUE(RESUME_PRODUCTION) && !isTRUE(FORCE_FRESH_PRODUCTION) &&
      file.exists(checkpoint_file)) {
    prior_checkpoint <- tryCatch(
      readRDS(checkpoint_file),
      error = function(e) NULL
    )
    was_complete_before <- is.list(prior_checkpoint) &&
      identical(
        prior_checkpoint$settings_signature,
        PRODUCTION_SETTINGS_SIGNATURE
      ) &&
      isTRUE(prior_checkpoint$complete)
  }
  
  production_results[[replicate_number]] <-
    run_one_production_replication(replicate_number)
  
  if (!isTRUE(production_results[[replicate_number]]$complete)) {
    stop("Production replication ", replicate_number, " did not complete.")
  }
  
  if (was_complete_before) {
    n_loaded_complete <- n_loaded_complete + 1L
  } else {
    n_computed_now <- n_computed_now + 1L
  }
  
  elapsed_sec <- as.numeric(difftime(
    Sys.time(),
    production_start,
    units = "secs"
  ))
  projected_total_sec <- elapsed_sec / replicate_number * N_REP
  
  cat(sprintf(
    paste0(
      "[%s] Completed PoE SAR replication %03d/%03d | ",
      "elapsed %s | projected total %s\n"
    ),
    timestamp_now(),
    replicate_number,
    N_REP,
    format_duration(elapsed_sec),
    format_duration(projected_total_sec)
  ))
}

production_end <- Sys.time()

## ============================================================================
## 17. COMBINE PRODUCTION RESULTS
## ============================================================================

if (any(!vapply(
  production_results,
  function(x) is.list(x) && isTRUE(x$complete),
  logical(1)
))) {
  stop("Not all 100 PoE SAR production replications completed.")
}

metrics_by_replicate <- do.call(
  rbind,
  lapply(production_results, `[[`, "metrics")
)
row.names(metrics_by_replicate) <- NULL

scores_all_nodes <- do.call(
  rbind,
  lapply(production_results, `[[`, "node_scores")
)
row.names(scores_all_nodes) <- NULL

geometry_summary <- do.call(
  rbind,
  lapply(production_results, `[[`, "geometry_summary")
)
row.names(geometry_summary) <- NULL

modularity_summary <- do.call(
  rbind,
  lapply(production_results, `[[`, "modularity_summary")
)
row.names(modularity_summary) <- NULL

prior_calibration_summary <- do.call(
  rbind,
  lapply(production_results, `[[`, "prior_calibration_summary")
)
row.names(prior_calibration_summary) <- NULL

sar_selection_summary <- do.call(
  rbind,
  lapply(production_results, `[[`, "sar_selection_summary")
)
row.names(sar_selection_summary) <- NULL

sar_selection_diagnostics <- do.call(
  rbind,
  lapply(production_results, `[[`, "sar_selection_diagnostics")
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

failure_parts <- Filter(
  Negate(is.null),
  lapply(production_results, `[[`, "failure")
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
  do.call(rbind, failure_parts)
}

if (nrow(metrics_by_replicate) != N_REP) {
  stop("The metric table does not contain exactly 100 rows.")
}
if (nrow(scores_all_nodes) != N_REP * N_NODES) {
  stop("The node-score table does not contain exactly 100 x 400 rows.")
}
if (nrow(prior_calibration_summary) != 2L * N_REP) {
  stop("The prior-calibration table does not contain exactly 2 x 100 rows.")
}
if (
  nrow(sar_selection_summary) != N_REP ||
  nrow(sar_selection_diagnostics) != N_REP * N_TRUE_POS
) {
  stop("The combined SAR-selection diagnostics have invalid dimensions.")
}

metrics_by_replicate <- metrics_by_replicate[
  order(metrics_by_replicate$replicate),
  ,
  drop = FALSE
]
scores_all_nodes <- scores_all_nodes[
  order(scores_all_nodes$replicate, scores_all_nodes$node),
  ,
  drop = FALSE
]
geometry_summary <- geometry_summary[
  order(geometry_summary$replicate),
  ,
  drop = FALSE
]
modularity_summary <- modularity_summary[
  order(modularity_summary$replicate),
  ,
  drop = FALSE
]
sar_selection_summary <- sar_selection_summary[
  order(sar_selection_summary$replicate),
  ,
  drop = FALSE
]
sar_selection_diagnostics <- sar_selection_diagnostics[
  order(
    sar_selection_diagnostics$replicate,
    sar_selection_diagnostics$node
  ),
  ,
  drop = FALSE
]
prior_calibration_summary <- prior_calibration_summary[
  order(
    prior_calibration_summary$replicate,
    match(prior_calibration_summary$expert, c("cov", "network"))
  ),
  ,
  drop = FALSE
]

write.csv(
  metrics_by_replicate,
  file.path(OUT_DIR, "PoE_metrics_by_replicate.csv"),
  row.names = FALSE
)
write.csv(
  scores_all_nodes,
  file.path(OUT_DIR, "PoE_all_node_probabilities.csv"),
  row.names = FALSE
)
write.csv(
  geometry_summary,
  file.path(OUT_DIR, "PoE_geometry_summary.csv"),
  row.names = FALSE
)
write.csv(
  modularity_summary,
  file.path(OUT_DIR, "PoE_modularity_summary.csv"),
  row.names = FALSE
)
write.csv(
  prior_calibration_summary,
  file.path(OUT_DIR, "PoE_prior_calibration_by_replicate.csv"),
  row.names = FALSE
)

write.csv(
  sar_selection_summary,
  file.path(OUT_DIR, "SAR_selection_summary_by_replicate.csv"),
  row.names = FALSE
)
write.csv(
  sar_selection_diagnostics,
  file.path(OUT_DIR, "SAR_selection_diagnostics_all_positives.csv"),
  row.names = FALSE
)
write.csv(
  sar_selection_mc_summary,
  file.path(OUT_DIR, "SAR_selection_Monte_Carlo_summary.csv"),
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

prior_calibration_distribution_summary <- do.call(
  rbind,
  lapply(c("cov", "network"), function(expert_name) {
    expert_rows <- prior_calibration_summary[
      prior_calibration_summary$expert == expert_name,
      ,
      drop = FALSE
    ]
    
    do.call(
      rbind,
      lapply(prior_summary_variables, function(variable_name) {
        statistics <- summary_stats(expert_rows[[variable_name]])
        data.frame(
          expert = expert_name,
          parameter = variable_name,
          n_rep = nrow(expert_rows),
          mean = unname(statistics["mean"]),
          sd = unname(statistics["sd"]),
          mc_se = unname(statistics["mc_se"]),
          q025 = unname(statistics["q025"]),
          median = unname(statistics["median"]),
          q975 = unname(statistics["q975"]),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
row.names(prior_calibration_distribution_summary) <- NULL

write.csv(
  prior_calibration_distribution_summary,
  file.path(OUT_DIR, "PoE_prior_calibration_distribution_summary.csv"),
  row.names = FALSE
)
write.csv(
  failure_log,
  file.path(OUT_DIR, "PoE_failure_log.csv"),
  row.names = FALSE
)

## ============================================================================
## 18. MONTE CARLO SUMMARIES
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
  Overall_Coverage_80 = "higher_with_interval_score",
  Hidden_Coverage_80 = "higher_with_interval_score",
  Overall_Interval_Score_80 = "lower",
  Hidden_Interval_Score_80 = "lower",
  Overall_Coverage_95 = "higher_with_interval_score",
  Hidden_Coverage_95 = "higher_with_interval_score",
  Overall_Interval_Score_95 = "lower",
  Hidden_Interval_Score_95 = "lower"
)

metric_summary_long <- do.call(
  rbind,
  lapply(main_metric_names, function(metric_name) {
    statistics <- summary_stats(metrics_by_replicate[[metric_name]])
    data.frame(
      scenario = SCENARIO_NAME,
      method = METHOD_NAME,
      metric = metric_name,
      direction = unname(metric_direction[metric_name]),
      n_rep = N_REP,
      mean = unname(statistics["mean"]),
      sd = unname(statistics["sd"]),
      mc_se = unname(statistics["mc_se"]),
      q025 = unname(statistics["q025"]),
      median = unname(statistics["median"]),
      q975 = unname(statistics["q975"]),
      stringsAsFactors = FALSE
    )
  })
)
row.names(metric_summary_long) <- NULL

write.csv(
  metric_summary_long,
  file.path(OUT_DIR, "PoE_metric_summary_long.csv"),
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
  file.path(OUT_DIR, "PoE_FINAL_MAIN_TABLE_means.csv"),
  row.names = FALSE
)

format_mean_sd <- function(mean_value, sd_value, digits = 3L) {
  if (!is.finite(mean_value)) {
    return("")
  }
  paste0(
    formatC(mean_value, format = "f", digits = digits),
    " (",
    formatC(sd_value, format = "f", digits = digits),
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
  Overall_Recall_at_45 = "Overall Recall@45",
  Hidden_Recall_at_15 = "Hidden Recall@15",
  Overall_Coverage_80 = "Overall Coverage 80%",
  Hidden_Coverage_80 = "Hidden Coverage 80%",
  Overall_Interval_Score_80 = "Overall Interval Score 80%",
  Hidden_Interval_Score_80 = "Hidden Interval Score 80%",
  Overall_Coverage_95 = "Overall Coverage 95%",
  Hidden_Coverage_95 = "Hidden Coverage 95%",
  Overall_Interval_Score_95 = "Overall Interval Score 95%",
  Hidden_Interval_Score_95 = "Hidden Interval Score 95%"
)

paper_table_mean_sd <- data.frame(
  Method = METHOD_NAME,
  stringsAsFactors = FALSE
)

for (metric_name in main_metric_names) {
  values <- metrics_by_replicate[[metric_name]]
  paper_table_mean_sd[[paper_column_labels[metric_name]]] <-
    format_mean_sd(
      mean(values, na.rm = TRUE),
      sd(values, na.rm = TRUE),
      digits = 3L
    )
}

write.csv(
  paper_table_mean_sd,
  file.path(OUT_DIR, "PoE_FINAL_PAPER_TABLE_mean_SD.csv"),
  row.names = FALSE
)

additional_metric_names <- c(
  "Hidden_at_15",
  "Overall_Interval_Width_80",
  "Hidden_Interval_Width_80",
  "Overall_Interval_Width_95",
  "Hidden_Interval_Width_95",
  "anchor_max_median",
  "anchor_Pr_over_090",
  "anchor_Pr_over_095",
  "anchor_Pr_over_099",
  "anchor_q99_median",
  "anchor_top5_median",
  "beta0_posterior_mean",
  "tau_cov_posterior_mean",
  "tau_network_posterior_mean",
  "ell_cov_posterior_mean",
  "sigma_cov_posterior_mean",
  "ell_network_posterior_mean",
  "sigma_network_posterior_mean",
  "network_mean_feature_count",
  "acc_blocked_f",
  "acc_blocked_cov",
  "acc_blocked_network",
  "acc_joint_tune",
  "acc_joint_sampling",
  "acc_f_refresh_sampling",
  "delta_final",
  "kappa_cov_final",
  "kappa_network_final",
  "joint_theta_scale_final",
  "runtime_sec"
)

additional_diagnostic_summary <- do.call(
  rbind,
  lapply(additional_metric_names, function(metric_name) {
    statistics <- summary_stats(metrics_by_replicate[[metric_name]])
    data.frame(
      method = METHOD_NAME,
      metric = metric_name,
      n_rep = N_REP,
      mean = unname(statistics["mean"]),
      sd = unname(statistics["sd"]),
      mc_se = unname(statistics["mc_se"]),
      q025 = unname(statistics["q025"]),
      median = unname(statistics["median"]),
      q975 = unname(statistics["q975"]),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  additional_diagnostic_summary,
  file.path(OUT_DIR, "PoE_additional_diagnostic_summary.csv"),
  row.names = FALSE
)

## ============================================================================
## 19. SAVE THE COMPLETE OUTPUT BUNDLE
## ============================================================================

total_production_runtime_sec <- as.numeric(difftime(
  production_end,
  production_start,
  units = "secs"
))

run_config <- list(
  script_version = SCRIPT_VERSION,
  evaluation_version = EVALUATION_VERSION,
  created_at = timestamp_now(),
  scenario = SCENARIO_NAME,
  method = METHOD_NAME,
  model_definition = paste0(
    "Two-expert ordinary Bernoulli-logistic PoE: covariate GP + one ",
    "network GP, with covariate and modularity mean terms."
  ),
  label_mechanism = list(
    type = "covariate-driven SAR",
    lambda = LAMBDA_SAR,
    z_clip = SAR_Z_CLIP,
    n_hidden_per_replication = N_HIDDEN,
    n_observed_positive_per_replication = N_OBSERVED_POS,
    selection = paste0(
      "Exactly 15 of the 45 true positives are sampled without replacement ",
      "with weights proportional to exp(-lambda*z_i), where z_i is the ",
      "oracle covariate fraud-likeness score standardized among the true ",
      "positives and clipped to [-3,3]."
    ),
    oracle_score_used_for_fitting = FALSE
  ),
  prior_calibration_mode = "exact_per_geometry",
  prior_calibration_rules = PRIOR_CALIBRATION_RULES,
  beta0_prior = beta0_prior_for_sampler,
  fixed_tau_cov_scale = TAU_SCALE_COV_ACTIVE,
  active_prior_rounding = ROUND_ACTIVE_PRIORS,
  active_prior_round_digits = PRIOR_ROUND_DIGITS,
  global_calibration_file = GLOBAL_CALIBRATION_FILE,
  global_calibration_hash = GLOBAL_CALIBRATION_HASH,
  n_pilot = N_PILOT,
  pilot_burnin = PILOT_BURNIN,
  pilot_tail_length = PILOT_TAIL_LENGTH,
  n_rep = N_REP,
  production_burnin = PRODUCTION_BURNIN,
  production_sampling = PRODUCTION_SAMPLING,
  production_thin = PRODUCTION_THIN,
  true_labels = TRUE_T,
  observed_labels =
    "replication-specific SAR labels stored in production_results$label_design",
  production_seed_table = production_seed_table,
  pilot_seed_table = pilot_seed_table,
  use_saved_competing_geometries = USE_SAVED_COMPETING_GEOMETRIES,
  sar_selection_summary = sar_selection_summary,
  sar_selection_Monte_Carlo_summary = sar_selection_mc_summary,
  n_loaded_complete = n_loaded_complete,
  n_computed_now = n_computed_now,
  total_production_runtime_sec = total_production_runtime_sec,
  predictive_interval_levels = PREDICTIVE_INTERVAL_LEVELS,
  predictive_label_construction = paste0(
    "For every retained posterior state and every node, draw ",
    "T_pred from Bernoulli(plogis(f)). Use central discrete empirical ",
    "quantiles with type=1."
  ),
  evaluation_truth = paste0(
    "Intervals and point metrics are compared with fixed TRUE_T: 0 for true ",
    "zeros and 1 for both observed and SAR-hidden positives. The ground truth ",
    "is never regenerated."
  ),
  H_intervals_used = FALSE,
  poisson_binomial_used = FALSE
)

saveRDS(
  run_config,
  file.path(OUT_DIR, "PoE_run_config.rds")
)

save(
  global_pilot_calibration,
  production_results,
  metrics_by_replicate,
  scores_all_nodes,
  geometry_summary,
  modularity_summary,
  sar_selection_summary,
  sar_selection_diagnostics,
  sar_selection_mc_summary,
  prior_calibration_summary,
  prior_calibration_distribution_summary,
  failure_log,
  metric_summary_long,
  main_table_means,
  paper_table_mean_sd,
  additional_diagnostic_summary,
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
  PRIOR_CALIBRATION_RULES,
  TAU_SCALE_COV_ACTIVE,
  run_config,
  PILOT_SETTINGS_SIGNATURE,
  PRODUCTION_SETTINGS_SIGNATURE,
  file = file.path(
    OUT_DIR,
    "SIMULATION_100_FULL_POE_SAR_SELF_CALIBRATING_COMPLETE_RESULTS.RData"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo_full_PoE_SAR_simulation.txt")
)

## ============================================================================
## 20. FINAL CONSOLE OUTPUT
## ============================================================================

cat("\n================ RUN ACCOUNTING: POE SAR ================\n")
cat("Pilot datasets          :", N_PILOT, "\n")
cat("Production replications :", N_REP, "\n")
cat("Computed in this run    :", n_computed_now, "\n")
cat("Loaded from checkpoint  :", n_loaded_complete, "\n")
cat("Production runtime      :",
    format_duration(total_production_runtime_sec), "\n")
cat("Prior calibration mode  : exact per pilot/replication\n")
cat("Application dependency  : none\n")

cat("\n================ SAR SELECTION MONTE CARLO SUMMARY ================\n")
print(sar_selection_mc_summary, row.names = FALSE)

cat("\n================ FINAL MONTE CARLO MEANS ================\n")
print(main_table_means, row.names = FALSE)

cat("\n================ FINAL PAPER TABLE: ALL METRICS ================\n")
print(paper_table_mean_sd, row.names = FALSE)


cat("\nAll outputs were written to:\n")
cat(normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE), "\n")
cat("\n================ DONE ================\n")
