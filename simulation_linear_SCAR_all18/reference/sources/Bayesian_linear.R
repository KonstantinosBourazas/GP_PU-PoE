## ============================================================================
## SIMULATION STUDY: FOUR BAYESIAN LINEAR MODELS ON LINEAR-GAUSSIAN + SBM SCAR DATA
## ----------------------------------------------------------------------------
## Models, in this order:
##   1) Bayesian Linear Cov
##   2) Bayesian Linear Cov + PU
##   3) Bayesian Linear Cov + Net
##   4) Bayesian Linear Cov + Net + PU
##
## The Bayesian models, priors, Polya-Gamma Gibbs updates, PU mechanism,
## five-pilot initialization, production MCMC, metrics and output tables are
## unchanged from the nonlinear-SCAR Bayesian-linear script.
##
## Only the data geometry is replaced by the calibrated linear SCAR design used
## by the competing-method study:
##
##   True zeros:
##       1:355
##
##   Hidden positives, observed as Y=0:
##       356:360  covariate signal only
##       361:365  network signal only
##       366:370  both covariate and network signal
##
##   Observed positives, Y=1:
##       371:380  covariate signal only
##       381:390  network signal only
##       391:400  both covariate and network signal
##
## Hence the original mixed-signal allocation is preserved exactly:
##   - hidden:   5 covariate-only + 5 network-only + 5 both;
##   - observed:10 covariate-only + 10 network-only + 10 both.
##
## Linear Gaussian covariates
## --------------------------
## The 30 covariate-signal nodes, namely covariate-only and both-signal
## positives, have ten independent N(0.5,1) variables. All remaining 370 nodes,
## including the 15 network-only positives, have ten independent N(0,1)
## variables.
##
## Bernoulli SBM
## -------------
## The 30 network-signal nodes, namely network-only and both-signal positives,
## form block 1. The remaining 370 nodes, including the 15 covariate-only
## positives, form block 0:
##
##   p_00 = 0.10,  p_01 = 0.08,  p_11 = 0.25.
##
## No Chinese-restaurant layer and no triadic closure are used.
##
## Five-pilot global initialization
## --------------------------------
## Five new linear-Gaussian/SBM pilot data sets are used. Old nonlinear pilot
## calibrations are not reused because the pilot data geometry has changed.
##
## Production geometry
## -------------------
## By default, each production replication loads the corresponding geometry
## saved by the competing-method linear-SCAR script. If the file is
## absent, the same geometry is regenerated from the identical deterministic
## covariate, network and modularity seeds.
##
## Production stage
## ----------------
##   - 100 independent data sets;
##   - 15,000 burn-in iterations;
##   - 40,000 posterior iterations;
##   - thinning 20, hence 2,000 retained draws per model and replication.
##
## The known latent truth TRUE_T is used for all evaluation metrics. Hidden
## positives therefore count as T=1 although their observed label is Y=0.
## ============================================================================

options(stringsAsFactors = FALSE)

## ============================================================================
## 0. PACKAGES
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

PG_BACKEND <- if (requireNamespace("BayesLogit", quietly = TRUE)) {
  "BayesLogit"
} else if (requireNamespace("pgdraw", quietly = TRUE)) {
  "pgdraw"
} else {
  stop(
    "A Polya-Gamma package is required. Install 'BayesLogit' ",
    "(preferred) or 'pgdraw'."
  )
}

suppressPackageStartupMessages({
  library(Matrix)
  library(igraph)
})

## ============================================================================
## 1. USER SETTINGS
## ============================================================================

SCRIPT_VERSION <- "bayesian_linear_linear_gaussian_SBM_SCAR_mixed_signal_5pilots_all_metrics_v1_1"

## The DGP has changed, so nonlinear pilot checkpoints must not be reused.
PILOT_CALIBRATION_VERSION <- paste0(
  "bayesian_linear_linear_gaussian_SBM_SCAR_",
  "mixed_signal_5pilots_v1_0"
)

MASTER_SEED <- 20260718L
N_PILOT <- 5L
N_REP <- 100L

PILOT_BURNIN <- 40000L
PILOT_TAIL_LENGTH <- 5000L

PRODUCTION_BURNIN <- 15000L
PRODUCTION_SAMPLING <- 40000L
PRODUCTION_THIN <- 20L

NU_TAU <- 3
INITIAL_ETA <- 0.20
ETA_EPS <- 1e-10

PG_OMEGA_FLOOR <- 1e-12
CHOLESKY_JITTER <- 1e-10
CHOLESKY_MAX_TRIES <- 8L

PILOT_PROGRESS_EVERY <- 5000L
PRODUCTION_BURNIN_PROGRESS_EVERY <- 5000L
PRODUCTION_SAMPLE_PROGRESS_EVERY <- 5000L

MOD_EIG_TOL <- 1e-8
MOD_MAX_Q <- Inf
MOD_K_MAX <- 20L

HIDDEN_TOP_K <- 15L
OVERALL_RECALL_K <- 45L

PREDICTIVE_INTERVAL_LEVELS <- c(0.80, 0.95)
PREDICTIVE_QUANTILE_TYPE <- 1L

PROBABILITY_EPS <- 1e-8

OUT_DIR <- file.path(
  getwd(),
  "SIMULATION_100_BAYESIAN_LINEAR_LINEAR_GAUSSIAN_SBM_SCAR_MIXED_SIGNAL_5_PILOTS"
)

PILOT_CHECKPOINT_DIR <- file.path(
  OUT_DIR,
  "pilot_checkpoints"
)

PRODUCTION_CHECKPOINT_DIR <- file.path(
  OUT_DIR,
  "production_checkpoints"
)

PILOT_GEOMETRY_DIR <- file.path(
  OUT_DIR,
  "pilot_geometries"
)

for (path_now in c(
  OUT_DIR,
  PILOT_CHECKPOINT_DIR,
  PRODUCTION_CHECKPOINT_DIR,
  PILOT_GEOMETRY_DIR
)) {
  dir.create(path_now, showWarnings = FALSE, recursive = TRUE)
}

GLOBAL_CALIBRATION_FILE <- file.path(
  OUT_DIR,
  "global_pilot_calibration_5_datasets.rds"
)

GLOBAL_CALIBRATION_RDATA <- file.path(
  OUT_DIR,
  "global_pilot_calibration_5_datasets.RData"
)

RESUME_PILOTS <- TRUE
RESUME_PRODUCTION <- TRUE
FORCE_FRESH_PILOTS <- FALSE
FORCE_FRESH_PRODUCTION <- FALSE
FAIL_IF_ANY_MODEL_FAILS <- TRUE
SAVE_PILOT_GEOMETRIES <- TRUE

## By default, use the geometries generated by the linear-SCAR
## competing-method script. This guarantees identical X, A and modularity data.
USE_SAVED_COMPETING_GEOMETRIES <- TRUE
COMPETING_OUT_DIR <- file.path(
  getwd(),
  "SIMULATION_100_COMPETING_METHODS_LINEAR_GAUSSIAN_SBM_SCAR_MIXED_SIGNAL"
)
COMPETING_GEOMETRY_DIR <- file.path(
  COMPETING_OUT_DIR,
  "generated_geometries"
)
REQUIRE_MATCHING_LINEAR_SBM_GEOMETRY <- TRUE

## Priors are constructed internally, with the same rules as the GP
## scripts. No quantity is imported from the application.
ROUND_ACTIVE_PRIORS <- TRUE
PRIOR_ROUND_DIGITS <- 2L

## Total latent residual SD budget and its allocation. A one-group model
## receives the whole budget, as in GP-Cov. A two-group model splits it
## equally, as in GP-PoE and GP-PU-PoE.
TOTAL_SD_CENTER <- 2
COV_VARIANCE_SHARE <- 0.50
PAIR_VARIANCE_SHARES <- c(
  cov = COV_VARIANCE_SHARE,
  network = 1 - COV_VARIANCE_SHARE
)

## Structured-mean half-t calibration.
HALFT_DF <- 3
MEAN_RMS_95_MULT <- 4.16

## Fixed beta0 and eta prior rules.
BETA0_PREV_RANGE_95 <- c(0.01, 0.10)
ETA_ALPHA <- 2
ETA_Q95 <- 0.50

if (!identical(as.integer(NU_TAU), as.integer(HALFT_DF))) {
  stop("NU_TAU and HALFT_DF must be identical.")
}
if (!isTRUE(all.equal(sum(PAIR_VARIANCE_SHARES), 1)) ||
    any(PAIR_VARIANCE_SHARES <= 0)) {
  stop("The two-group variance allocation is inconsistent.")
}

## Production seeds match the competing-method script for
## covariates, networks and modularity features.
PRODUCTION_MCMC_SEED_BASE <- 810000000L
PRODUCTION_PREDICTIVE_SEED_BASE <- 910000000L
PILOT_MCMC_SEED_BASE <- 710000000L

MODEL_SEED_OFFSETS <- c(
  cov = 0L,
  cov_pu = 1000000L,
  cov_net = 2000000L,
  cov_net_pu = 3000000L
)

## ============================================================================
## 1A. FIXED DATA DESIGN
## ============================================================================

N_NODES <- 400L
TRUE_ZERO_IDX <- 1:355

HIDDEN_COV_ONLY_IDX <- 356:360
HIDDEN_NET_ONLY_IDX <- 361:365
HIDDEN_BOTH_IDX <- 366:370

OBSERVED_COV_ONLY_IDX <- 371:380
OBSERVED_NET_ONLY_IDX <- 381:390
OBSERVED_BOTH_IDX <- 391:400

HIDDEN_POS_IDX <- c(
  HIDDEN_COV_ONLY_IDX,
  HIDDEN_NET_ONLY_IDX,
  HIDDEN_BOTH_IDX
)

OBSERVED_POS_IDX <- c(
  OBSERVED_COV_ONLY_IDX,
  OBSERVED_NET_ONLY_IDX,
  OBSERVED_BOTH_IDX
)

TRUE_POS_IDX <- c(
  HIDDEN_POS_IDX,
  OBSERVED_POS_IDX
)

COVARIATE_SIGNAL_IDX <- c(
  HIDDEN_COV_ONLY_IDX,
  HIDDEN_BOTH_IDX,
  OBSERVED_COV_ONLY_IDX,
  OBSERVED_BOTH_IDX
)

NETWORK_SIGNAL_IDX <- c(
  HIDDEN_NET_ONLY_IDX,
  HIDDEN_BOTH_IDX,
  OBSERVED_NET_ONLY_IDX,
  OBSERVED_BOTH_IDX
)

N_TRUE_POS <- length(TRUE_POS_IDX)
N_HIDDEN <- length(HIDDEN_POS_IDX)
N_OBSERVED_POS <- length(OBSERVED_POS_IDX)
N_TRUE_ZERO <- length(TRUE_ZERO_IDX)
N_UNLABELED <- N_TRUE_ZERO + N_HIDDEN

N_COVARIATE_SIGNAL <- length(COVARIATE_SIGNAL_IDX)
N_NETWORK_SIGNAL <- length(NETWORK_SIGNAL_IDX)
N_BOTH_SIGNAL <- length(
  intersect(
    COVARIATE_SIGNAL_IDX,
    NETWORK_SIGNAL_IDX
  )
)
N_NETWORK_BACKGROUND <- N_NODES - N_NETWORK_SIGNAL

TRUE_T <- integer(N_NODES)
TRUE_T[TRUE_POS_IDX] <- 1L

OBSERVED_Y <- integer(N_NODES)
OBSERVED_Y[OBSERVED_POS_IDX] <- 1L

NETWORK_BLOCK <- integer(N_NODES)
NETWORK_BLOCK[NETWORK_SIGNAL_IDX] <- 1L

SIGNAL_TYPE <- rep("true zero", N_NODES)
SIGNAL_TYPE[HIDDEN_COV_ONLY_IDX] <- "hidden: covariates only"
SIGNAL_TYPE[HIDDEN_NET_ONLY_IDX] <- "hidden: network only"
SIGNAL_TYPE[HIDDEN_BOTH_IDX] <- "hidden: both"
SIGNAL_TYPE[OBSERVED_COV_ONLY_IDX] <- "observed: covariates only"
SIGNAL_TYPE[OBSERVED_NET_ONLY_IDX] <- "observed: network only"
SIGNAL_TYPE[OBSERVED_BOTH_IDX] <- "observed: both"

SCENARIO_NAME <- paste0(
  "Bayesian linear models on linear Gaussian covariates and Bernoulli SBM, ",
  "fixed mixed-signal SCAR design"
)

## ============================================================================
## 1B. LINEAR GAUSSIAN COVARIATE DGP SETTINGS
## ============================================================================

N_COVARIATES <- 10L
COVARIATE_BACKGROUND_MEAN <- 0.00
COVARIATE_SIGNAL_MEAN <- 0.50
COMMON_SD <- 1.00

## Backward-compatible aliases used only in saved metadata.
ZERO_MEAN <- COVARIATE_BACKGROUND_MEAN
POSITIVE_MEAN <- COVARIATE_SIGNAL_MEAN

## ============================================================================
## 1C. BERNOULLI SBM NETWORK DGP SETTINGS
## ============================================================================

P_00 <- 0.10
P_01 <- 0.08
P_11 <- 0.25


TRIAD_CLOSURE_PROB <- 0.00
TRIAD_MAX_NEW_EDGES <- Inf

## ============================================================================
## 1D. INTERNALLY CONSTRUCTED PRIORS
## ============================================================================
## The beta0 and eta rules are deterministic. The structured-mean half-t
## scales are calibrated from the geometry of the data, as in the GP
## scripts, so the network scale is computed per replication inside
## make_model_specs().

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

beta0_low <- logit(BETA0_PREV_RANGE_95[1L])
beta0_high <- logit(BETA0_PREV_RANGE_95[2L])

beta0_prior_for_sampler <- list(
  mean = round_prior(0.5 * (beta0_low + beta0_high)),
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

## Every covariate matrix has ten columns standardized with sample SD = 1.
## Therefore sum(X^2)/n = 10 * (n - 1)/n exactly up to floating-point error.
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

EFFECTIVE_SD_CENTER_SINGLE <- TOTAL_SD_CENTER
EFFECTIVE_SD_CENTER_PAIR <- TOTAL_SD_CENTER * sqrt(PAIR_VARIANCE_SHARES)

if (!isTRUE(all.equal(
  sqrt(sum(EFFECTIVE_SD_CENTER_PAIR^2)),
  TOTAL_SD_CENTER
))) {
  stop("The two-group effective SD allocation is inconsistent.")
}

half_t_scale_from_geometry <- function(effective_sd_center, b) {
  (MEAN_RMS_95_MULT * effective_sd_center) / (HALF_T_Q95_UNIT * b)
}

TAU_SCALE_COV_SINGLE_UNROUNDED <- half_t_scale_from_geometry(
  EFFECTIVE_SD_CENTER_SINGLE,
  B_COV_FIXED
)

TAU_SCALE_COV_SINGLE_ACTIVE <- if (isTRUE(ROUND_ACTIVE_PRIORS)) {
  round_prior(TAU_SCALE_COV_SINGLE_UNROUNDED)
} else {
  TAU_SCALE_COV_SINGLE_UNROUNDED
}

if (!is.finite(TAU_SCALE_COV_SINGLE_ACTIVE) ||
    TAU_SCALE_COV_SINGLE_ACTIVE <= 0) {
  stop("The one-group covariate half-t scale is invalid.")
}

TAU2_COV_SINGLE_DEFAULT <- TAU_SCALE_COV_SINGLE_ACTIVE^2

## Geometry-specific half-t scales for the two-group models. The covariate
## factor is fixed by standardization; the network factor depends on the
## modularity dimension and is therefore recomputed for every data set.
geometry_half_t_scales <- function(X_cov, Z_network) {
  X_cov <- as.matrix(X_cov)
  Z_network <- as.matrix(Z_network)

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

  tau_scale_unrounded <- c(
    cov = half_t_scale_from_geometry(
      unname(EFFECTIVE_SD_CENTER_PAIR["cov"]),
      B_COV_FIXED
    ),
    network = half_t_scale_from_geometry(
      unname(EFFECTIVE_SD_CENTER_PAIR["network"]),
      b_network
    )
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

  list(
    tau_scale = tau_scale_active,
    tau2_default = tau_scale_active^2,
    tau_scale_unrounded = tau_scale_unrounded,
    b_network = b_network
  )
}

## Recorded in the run signatures so that checkpoints produced under the
## previous application-based priors are invalidated.
TAU_SCALE_RULE <- list(
  source = "internal geometry calibration, no application quantities",
  total_sd_center = TOTAL_SD_CENTER,
  pair_variance_shares = PAIR_VARIANCE_SHARES,
  single_effective_sd_center = EFFECTIVE_SD_CENTER_SINGLE,
  halft_df = HALFT_DF,
  mean_rms_95_mult = MEAN_RMS_95_MULT,
  round_active_priors = ROUND_ACTIVE_PRIORS,
  prior_round_digits = PRIOR_ROUND_DIGITS,
  b_cov_fixed = B_COV_FIXED,
  tau_scale_cov_single = TAU_SCALE_COV_SINGLE_ACTIVE,
  network_scale = "recomputed per replication from the modularity dimension"
)


## ============================================================================
## 1E. REPRODUCIBLE SEED TABLES
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
  stringsAsFactors = FALSE
)

all_seed_values <- c(
  unlist(production_seed_table[-1L]),
  unlist(pilot_seed_table[-1L]),
  PRODUCTION_MCMC_SEED_BASE,
  PRODUCTION_PREDICTIVE_SEED_BASE,
  PILOT_MCMC_SEED_BASE,
  MODEL_SEED_OFFSETS
)

if (any(all_seed_values > .Machine$integer.max)) {
  stop("At least one deterministic seed exceeds the R integer range.")
}

write.csv(
  production_seed_table,
  file.path(OUT_DIR, "production_seed_table.csv"),
  row.names = FALSE
)

write.csv(
  pilot_seed_table,
  file.path(OUT_DIR, "pilot_seed_table.csv"),
  row.names = FALSE
)

## ============================================================================
## 2. DESIGN AND MCMC CHECKS
## ============================================================================

if (
  N_NODES != 400L ||
  N_TRUE_POS != 45L ||
  N_HIDDEN != 15L ||
  N_OBSERVED_POS != 30L ||
  N_TRUE_ZERO != 355L ||
  N_UNLABELED != 370L
) {
  stop("The fixed 355/15/30 simulation design has changed.")
}

if (
  !identical(TRUE_ZERO_IDX, 1:355) ||
  !identical(HIDDEN_COV_ONLY_IDX, 356:360) ||
  !identical(HIDDEN_NET_ONLY_IDX, 361:365) ||
  !identical(HIDDEN_BOTH_IDX, 366:370) ||
  !identical(OBSERVED_COV_ONLY_IDX, 371:380) ||
  !identical(OBSERVED_NET_ONLY_IDX, 381:390) ||
  !identical(OBSERVED_BOTH_IDX, 391:400) ||
  !identical(HIDDEN_POS_IDX, 356:370) ||
  !identical(OBSERVED_POS_IDX, 371:400) ||
  !identical(TRUE_POS_IDX, 356:400)
) {
  stop("The fixed mixed-signal row-index design is inconsistent.")
}

if (
  length(HIDDEN_COV_ONLY_IDX) != 5L ||
  length(HIDDEN_NET_ONLY_IDX) != 5L ||
  length(HIDDEN_BOTH_IDX) != 5L ||
  length(OBSERVED_COV_ONLY_IDX) != 10L ||
  length(OBSERVED_NET_ONLY_IDX) != 10L ||
  length(OBSERVED_BOTH_IDX) != 10L
) {
  stop("The required hidden 5+5+5 and observed 10+10+10 split has changed.")
}

if (
  N_COVARIATE_SIGNAL != 30L ||
  N_NETWORK_SIGNAL != 30L ||
  N_BOTH_SIGNAL != 15L
) {
  stop("The design requires 30 covariate, 30 network and 15 both-signal nodes.")
}

if (
  !identical(
    intersect(COVARIATE_SIGNAL_IDX, NETWORK_SIGNAL_IDX),
    c(HIDDEN_BOTH_IDX, OBSERVED_BOTH_IDX)
  )
) {
  stop("The both-signal allocation is inconsistent.")
}

if (
  !identical(which(TRUE_T == 1L), TRUE_POS_IDX) ||
  !identical(which(OBSERVED_Y == 1L), OBSERVED_POS_IDX) ||
  !identical(
    which(TRUE_T == 1L & OBSERVED_Y == 0L),
    HIDDEN_POS_IDX
  )
) {
  stop("TRUE_T and OBSERVED_Y do not match the requested SCAR design.")
}

if (!identical(which(NETWORK_BLOCK == 1L), NETWORK_SIGNAL_IDX)) {
  stop("NETWORK_BLOCK does not match NETWORK_SIGNAL_IDX.")
}

if (
  N_COVARIATES != 10L ||
  COVARIATE_BACKGROUND_MEAN != 0 ||
  COVARIATE_SIGNAL_MEAN != 0.5 ||
  COMMON_SD != 1
) {
  stop("The requested covariate DGP is ten N(0,1)/N(0.5,1) variables.")
}


EXPECTED_DEGREE_NETWORK_BACKGROUND <-
  (N_NETWORK_BACKGROUND - 1L) * P_00 +
  N_NETWORK_SIGNAL * P_01

EXPECTED_DEGREE_NETWORK_SIGNAL <-
  N_NETWORK_BACKGROUND * P_01 +
  (N_NETWORK_SIGNAL - 1L) * P_11

EXPECTED_EDGE_COUNTS_NETWORK_BLOCKS <- c(
  "00" = choose(N_NETWORK_BACKGROUND, 2L) * P_00,
  "01" = N_NETWORK_BACKGROUND * N_NETWORK_SIGNAL * P_01,
  "11" = choose(N_NETWORK_SIGNAL, 2L) * P_11
)

if (N_PILOT != 5L) {
  stop("This script requires exactly five independent pilots.")
}

if (
  PILOT_BURNIN <= 0L ||
  PILOT_TAIL_LENGTH <= 0L ||
  PILOT_TAIL_LENGTH > PILOT_BURNIN
) {
  stop("Invalid pilot warm-up or tail length.")
}

if (
  PRODUCTION_BURNIN < 0L ||
  PRODUCTION_SAMPLING <= 0L ||
  PRODUCTION_THIN <= 0L ||
  PRODUCTION_SAMPLING %% PRODUCTION_THIN != 0L
) {
  stop("Invalid production MCMC lengths.")
}

if (
  length(PREDICTIVE_INTERVAL_LEVELS) != 2L ||
  !isTRUE(all.equal(PREDICTIVE_INTERVAL_LEVELS, c(0.80, 0.95))) ||
  any(!is.finite(PREDICTIVE_INTERVAL_LEVELS)) ||
  any(PREDICTIVE_INTERVAL_LEVELS <= 0) ||
  any(PREDICTIVE_INTERVAL_LEVELS >= 1)
) {
  stop("PREDICTIVE_INTERVAL_LEVELS must equal c(0.80,0.95).")
}

if (!identical(as.integer(PREDICTIVE_QUANTILE_TYPE), 1L)) {
  stop("PREDICTIVE_QUANTILE_TYPE must equal 1L for binary intervals.")
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

cat("\n================ BAYESIAN LINEAR-SCAR SIMULATION DESIGN ================\n")
cat("Prior construction        : internal geometry calibration\n")
cat("Polya-Gamma backend           :", PG_BACKEND, "\n")
cat("True zeros / positives        :", N_TRUE_ZERO, "/", N_TRUE_POS, "\n")
cat("Observed / hidden positives   :", N_OBSERVED_POS, "/", N_HIDDEN, "\n")
cat(
  "Hidden cov/net/both          :",
  length(HIDDEN_COV_ONLY_IDX), "/",
  length(HIDDEN_NET_ONLY_IDX), "/",
  length(HIDDEN_BOTH_IDX), "\n"
)
cat(
  "Observed cov/net/both        :",
  length(OBSERVED_COV_ONLY_IDX), "/",
  length(OBSERVED_NET_ONLY_IDX), "/",
  length(OBSERVED_BOTH_IDX), "\n"
)
cat(
  "Covariate/network/both signal:",
  N_COVARIATE_SIGNAL, "/",
  N_NETWORK_SIGNAL, "/",
  N_BOTH_SIGNAL, "\n"
)
cat(
  "Covariate background/signal  : N(0,1) / N(0.5,1),",
  N_COVARIATES, "variables\n"
)
cat("SBM p00/p01/p11              :", P_00, "/", P_01, "/", P_11, "\n")
cat("Network block sizes 0/1      :", N_NETWORK_BACKGROUND, "/", N_NETWORK_SIGNAL, "\n")
cat(
  "Expected degree block 0/1    :",
  round(EXPECTED_DEGREE_NETWORK_BACKGROUND, 3), "/",
  round(EXPECTED_DEGREE_NETWORK_SIGNAL, 3), "\n"
)
cat("Triadic closure              :", TRIAD_CLOSURE_PROB, "\n")
cat("Pilot datasets               :", N_PILOT, "\n")
cat("Pilot warm-up                :", PILOT_BURNIN, "\n")
cat("Pilot pooled tail            :", PILOT_TAIL_LENGTH, "per pilot\n")
cat("Production replications      :", N_REP, "\n")
cat("Production burn-in           :", PRODUCTION_BURNIN, "\n")
cat("Production sampling          :", PRODUCTION_SAMPLING, "\n")
cat("Production thinning          :", PRODUCTION_THIN, "\n")
cat("Retained draws/run           :", PRODUCTION_SAMPLING / PRODUCTION_THIN, "\n")
cat(
  "Predictive interval levels   :",
  paste(PREDICTIVE_INTERVAL_LEVELS, collapse = ", "), "\n"
)
cat("Predictive quantile type     :", PREDICTIVE_QUANTILE_TYPE, "\n")
cat("beta0 prior (mean, sd)    :", beta0_prior_for_sampler$mean,
    ",", beta0_prior_for_sampler$sd, "\n")
cat("eta prior Beta(a, b)      :", eta_prior_for_sampler$a_eta,
    ",", eta_prior_for_sampler$b_eta, "\n")
cat("tau half-t df             :", HALFT_DF, "\n")
cat("Total latent SD budget    :", TOTAL_SD_CENTER, "\n")
cat("One-group tau_cov scale   :", TAU_SCALE_COV_SINGLE_ACTIVE, "\n")
cat("Two-group effective SDs   :\n")
print(EFFECTIVE_SD_CENTER_PAIR)
cat("Network tau scale         : recomputed per replication\n")

## ============================================================================
## 3. GENERAL HELPERS
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

sanitize_network_adjacency <- function(A) {
  A <- as.matrix(A)
  A[!is.finite(A)] <- 0
  A <- 1L * ((A + t(A)) > 0)
  diag(A) <- 0L
  A
}

chol_safe <- function(
    K,
    jitter = CHOLESKY_JITTER,
    max_tries = CHOLESKY_MAX_TRIES
) {
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

draw_inverse_gamma_one <- function(shape, scale, max_tries = 20L) {
  if (!is.finite(shape) || shape <= 0 ||
      !is.finite(scale) || scale <= 0) {
    stop("Invalid inverse-gamma parameters.")
  }
  for (ii in seq_len(as.integer(max_tries))) {
    g <- rgamma(1L, shape = shape, rate = scale)
    out <- 1 / g
    if (is.finite(out) && out > 0) return(out)
  }
  stop("Could not draw a finite inverse-gamma value.")
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

validate_groups <- function(groups, p_features) {
  if (!is.list(groups) || length(groups) < 1L || is.null(names(groups))) {
    stop("groups must be a named non-empty list.")
  }
  idx_raw <- as.integer(unlist(groups, use.names = FALSE))
  if (length(idx_raw) != p_features ||
      anyDuplicated(idx_raw) ||
      !identical(sort(idx_raw), seq_len(p_features))) {
    stop("The coefficient groups must partition all feature columns exactly.")
  }
  invisible(TRUE)
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

validate_probability_vector <- function(probabilities, n_expected, method_label) {
  probabilities <- as.numeric(probabilities)
  if (length(probabilities) != n_expected) {
    stop(method_label, " returned the wrong number of probabilities.")
  }
  if (any(!is.finite(probabilities))) {
    stop(method_label, " returned non-finite probabilities.")
  }
  tolerance <- 1e-8
  if (any(probabilities < -tolerance) ||
      any(probabilities > 1 + tolerance)) {
    stop(method_label, " returned probabilities outside [0,1].")
  }
  pmin(pmax(probabilities, 0), 1)
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

hidden_logloss <- function(
    probabilities,
    hidden_idx = HIDDEN_POS_IDX,
    eps = PROBABILITY_EPS
) {
  probabilities <- pmin(
    pmax(as.numeric(probabilities[hidden_idx]), eps),
    1 - eps
  )
  -mean(log(probabilities))
}

hidden_brier <- function(
    probabilities,
    hidden_idx = HIDDEN_POS_IDX
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
    hidden_idx = HIDDEN_POS_IDX,
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
  
  c(
    mean = mean(x),
    sd = stats::sd(x),
    mc_se = stats::sd(x) / sqrt(length(x)),
    q025 = unname(stats::quantile(x, 0.025, type = 7)),
    median = stats::median(x),
    q975 = unname(stats::quantile(x, 0.975, type = 7))
  )
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

file_md5_or_na <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path)[1L])
}


## ============================================================================
## 4. LINEAR GAUSSIAN COVARIATE DATA-GENERATING PROCESS
## ============================================================================

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
  
  if (
    !identical(
      signal_idx,
      COVARIATE_SIGNAL_IDX
    )
  ) {
    stop(
      "The Gaussian signal set must equal COVARIATE_SIGNAL_IDX."
    )
  }
  
  if (length(signal_idx) != 30L) {
    stop(
      "The covariate DGP requires exactly 30 signal nodes."
    )
  }
  
  background_idx <- setdiff(
    seq_len(N_NODES),
    signal_idx
  )
  
  expected_background_idx <- sort(
    c(
      TRUE_ZERO_IDX,
      HIDDEN_NET_ONLY_IDX,
      OBSERVED_NET_ONLY_IDX
    )
  )
  
  if (
    !identical(
      background_idx,
      expected_background_idx
    )
  ) {
    stop(
      "The covariate-background set must contain true zeros and network-only positives."
    )
  }
  
  X <- matrix(
    NA_real_,
    nrow = N_NODES,
    ncol = N_COVARIATES,
    dimnames = list(
      NULL,
      paste0(
        "X",
        seq_len(N_COVARIATES)
      )
    )
  )
  
  X[background_idx, ] <- matrix(
    stats::rnorm(
      length(background_idx) *
        N_COVARIATES,
      mean = COVARIATE_BACKGROUND_MEAN,
      sd = COMMON_SD
    ),
    nrow = length(background_idx),
    ncol = N_COVARIATES
  )
  
  X[signal_idx, ] <- matrix(
    stats::rnorm(
      length(signal_idx) *
        N_COVARIATES,
      mean = COVARIATE_SIGNAL_MEAN,
      sd = COMMON_SD
    ),
    nrow = length(signal_idx),
    ncol = N_COVARIATES
  )
  
  if (
    any(!is.finite(X)) ||
    !all(
      dim(X) ==
      c(
        N_NODES,
        N_COVARIATES
      )
    )
  ) {
    stop(
      "The generated Gaussian covariate matrix is invalid."
    )
  }
  
  standardized <- standardize_columns(
    X
  )
  
  covariate_profile <- rep(
    "covariate background: N(0,1)^10",
    N_NODES
  )
  
  covariate_profile[signal_idx] <-
    "covariate signal: N(0.5,1)^10"
  
  list(
    seed = seed,
    X_raw = X,
    X_std = standardized$scaled,
    center = standardized$center,
    scale = standardized$scale,
    signal_idx = signal_idx,
    background_idx = background_idx,
    profile_1_idx = integer(0),
    profile_2_idx = integer(0),
    covariate_profile = covariate_profile
  )
}

## ============================================================================
## 5. TWO-BLOCK BERNOULLI SBM DATA-GENERATING PROCESS
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
    nrow(
      as.matrix(pairs)
    ) == 0L
  ) {
    return(A)
  }
  
  pairs <- as.matrix(
    pairs
  )
  
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

sample_pair_set <- function(
    pairs,
    probability
) {
  pairs <- as.matrix(
    pairs
  )
  
  if (
    nrow(pairs) == 0L ||
    probability <= 0
  ) {
    return(
      matrix(
        integer(0),
        ncol = 2L
      )
    )
  }
  
  if (probability >= 1) {
    return(pairs)
  }
  
  keep <- stats::runif(
    nrow(pairs)
  ) <
    probability
  
  pairs[
    keep,
    ,
    drop = FALSE
  ]
}

build_bernoulli_sbm <- function(
    background_idx,
    signal_idx,
    p_00 = P_00,
    p_01 = P_01,
    p_11 = P_11
) {
  A <- matrix(
    0L,
    N_NODES,
    N_NODES
  )
  
  pairs_00 <- t(
    utils::combn(
      background_idx,
      2L
    )
  )
  
  pairs_01 <- cbind(
    rep(
      background_idx,
      each = length(signal_idx)
    ),
    rep(
      signal_idx,
      times = length(background_idx)
    )
  )
  
  pairs_11 <- t(
    utils::combn(
      signal_idx,
      2L
    )
  )
  
  selected_00 <- sample_pair_set(
    pairs_00,
    p_00
  )
  
  selected_01 <- sample_pair_set(
    pairs_01,
    p_01
  )
  
  selected_11 <- sample_pair_set(
    pairs_11,
    p_11
  )
  
  A <- add_edges(
    A,
    selected_00
  )
  
  A <- add_edges(
    A,
    selected_01
  )
  
  A <- add_edges(
    A,
    selected_11
  )
  
  A <- sanitize_network_adjacency(
    A
  )
  
  list(
    A = A,
    
    probabilities = c(
      "00" = p_00,
      "01" = p_01,
      "11" = p_11
    ),
    
    expected_edge_counts = c(
      "00" = nrow(pairs_00) * p_00,
      "01" = nrow(pairs_01) * p_01,
      "11" = nrow(pairs_11) * p_11
    ),
    
    realized_edge_counts = c(
      "00" = nrow(selected_00),
      "01" = nrow(selected_01),
      "11" = nrow(selected_11)
    ),
    
    edges = edge_count(A)
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
  
  if (closure_prob <= 0) {
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
    nrow(candidate_pairs) == 0L
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
      1 -
        closure_prob
    )^common_neighbor_count
  
  keep <- stats::runif(
    length(pair_probability)
  ) <
    pair_probability
  
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

edge_block_counts <- function(
    A,
    block_indicator = NETWORK_BLOCK
) {
  edge_index <- which(
    upper.tri(A) &
      A > 0,
    arr.ind = TRUE
  )
  
  if (
    nrow(edge_index) == 0L
  ) {
    return(
      c(
        "00" = 0L,
        "01" = 0L,
        "11" = 0L
      )
    )
  }
  
  endpoint_sum <-
    block_indicator[
      edge_index[, 1L]
    ] +
    block_indicator[
      edge_index[, 2L]
    ]
  
  c(
    "00" = sum(
      endpoint_sum == 0L
    ),
    "01" = sum(
      endpoint_sum == 1L
    ),
    "11" = sum(
      endpoint_sum == 2L
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
  
  if (
    !identical(
      signal_idx,
      NETWORK_SIGNAL_IDX
    )
  ) {
    stop(
      "The SBM signal block must equal NETWORK_SIGNAL_IDX."
    )
  }
  
  if (length(signal_idx) != 30L) {
    stop(
      "The SBM signal block must contain exactly 30 nodes."
    )
  }
  
  background_idx <- setdiff(
    seq_len(N_NODES),
    signal_idx
  )
  
  expected_background_idx <- sort(
    c(
      TRUE_ZERO_IDX,
      HIDDEN_COV_ONLY_IDX,
      OBSERVED_COV_ONLY_IDX
    )
  )
  
  if (
    !identical(
      background_idx,
      expected_background_idx
    )
  ) {
    stop(
      "The SBM background block must contain true zeros and covariate-only positives."
    )
  }
  
  sbm <- build_bernoulli_sbm(
    background_idx = background_idx,
    signal_idx = signal_idx,
    p_00 = P_00,
    p_01 = P_01,
    p_11 = P_11
  )
  
  A_preclosure <- sbm$A
  
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
  
  final_block_counts <- edge_block_counts(
    A_final,
    NETWORK_BLOCK
  )
  
  expected_degree_background <-
    (
      length(background_idx) -
        1L
    ) *
    P_00 +
    length(signal_idx) *
    P_01
  
  expected_degree_signal <-
    length(background_idx) *
    P_01 +
    (
      length(signal_idx) -
        1L
    ) *
    P_11
  
  summary <- data.frame(
    n_nodes = N_NODES,
    
    n_network_signal =
      length(signal_idx),
    
    n_network_background =
      length(background_idx),
    
    p_00 = P_00,
    p_01 = P_01,
    p_11 = P_11,
    
    expected_degree_network_background =
      expected_degree_background,
    
    expected_degree_network_signal =
      expected_degree_signal,
    
    expected_edges_network_00 =
      unname(
        sbm$expected_edge_counts["00"]
      ),
    
    expected_edges_network_01 =
      unname(
        sbm$expected_edge_counts["01"]
      ),
    
    expected_edges_network_11 =
      unname(
        sbm$expected_edge_counts["11"]
      ),
    
    realized_edges_network_00 =
      unname(
        final_block_counts["00"]
      ),
    
    realized_edges_network_01 =
      unname(
        final_block_counts["01"]
      ),
    
    realized_edges_network_11 =
      unname(
        final_block_counts["11"]
      ),
    
    edges_sbm = sbm$edges,
    edges_closure = closure$n_edges,
    edges_final = edge_count(A_final),
    
    closure_candidate_pairs =
      closure$n_candidates,
    
    mean_degree = mean(degree),
    sd_degree = stats::sd(degree),
    
    mean_degree_network_background =
      mean(
        degree[background_idx]
      ),
    
    sd_degree_network_background =
      stats::sd(
        degree[background_idx]
      ),
    
    mean_degree_network_signal =
      mean(
        degree[signal_idx]
      ),
    
    sd_degree_network_signal =
      stats::sd(
        degree[signal_idx]
      ),
    
    min_degree = min(degree),
    max_degree = max(degree),
    
    n_components =
      component_info$no,
    
    largest_component =
      max(
        component_info$csize
      ),
    
    stringsAsFactors = FALSE
  )
  
  list(
    seed = seed,
    A = A_final,
    A_sbm = sbm$A,
    A_closure = closure$A,
    signal_idx = signal_idx,
    background_idx = background_idx,
    block_indicator = NETWORK_BLOCK,
    degree = degree,
    summary = summary
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

## ============================================================================
## 7. FEATURE BUNDLES AND MODEL SPECIFICATIONS
## ============================================================================

make_generated_feature_bundle <- function(
    seed_row,
    source_label
) {
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
  
  X_cov <- as.matrix(covariate_set$X_std)
  Z_net <- as.matrix(modularity$Z)
  
  if (is.null(colnames(X_cov))) {
    colnames(X_cov) <- paste0("X", seq_len(ncol(X_cov)))
  }
  
  if (is.null(colnames(Z_net))) {
    colnames(Z_net) <- paste0("network_mod", seq_len(ncol(Z_net)))
  }
  
  X_cov_net <- cbind(X_cov, Z_net)
  colnames(X_cov_net) <- make.unique(colnames(X_cov_net))
  
  list(
    source = source_label,
    X_cov = X_cov,
    Z_net = Z_net,
    X_cov_net = X_cov_net,
    A_network = network_set$A,
    network_summary = network_set$summary,
    modularity_q = modularity$q,
    modularity_selected_eigenvalues =
      modularity$selected_eigenvalues,
    seeds = seed_row,
    covariate_profile =
      covariate_set$covariate_profile
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
    
    
    
    X_cov <- as.matrix(geometry$X_std)
    A_network <- sanitize_network_adjacency(
      geometry$A_network
    )
    
    if (
      !all(dim(X_cov) == c(N_NODES, N_COVARIATES)) ||
      !all(dim(A_network) == c(N_NODES, N_NODES))
    ) {
      stop("Saved competing geometry has incompatible dimensions.")
    }
    
    if (any(!is.finite(X_cov))) {
      stop("Saved competing covariates contain non-finite values.")
    }
    
    if (!is.null(geometry$Z_network_modularity)) {
      Z_net <- as.matrix(
        geometry$Z_network_modularity
      )
      
      modularity_q <- ncol(Z_net)
      
      modularity_selected_eigenvalues <-
        geometry$modularity_selected_eigenvalues
    } else {
      modularity <- modularity_features_single_network(
        A = A_network,
        seed = seed_row$modularity_seed,
        eig_tol = MOD_EIG_TOL,
        max_q = MOD_MAX_Q,
        k_max = MOD_K_MAX
      )
      
      Z_net <- as.matrix(modularity$Z)
      modularity_q <- modularity$q
      
      modularity_selected_eigenvalues <-
        modularity$selected_eigenvalues
    }
    
    if (
      nrow(Z_net) != N_NODES ||
      ncol(Z_net) < 1L ||
      any(!is.finite(Z_net))
    ) {
      stop("Saved or reconstructed network features are invalid.")
    }
    
    if (is.null(colnames(X_cov))) {
      colnames(X_cov) <- paste0(
        "X",
        seq_len(ncol(X_cov))
      )
    }
    
    if (is.null(colnames(Z_net))) {
      colnames(Z_net) <- paste0(
        "network_mod",
        seq_len(ncol(Z_net))
      )
    }
    
    X_cov_net <- cbind(
      X_cov,
      Z_net
    )
    
    colnames(X_cov_net) <- make.unique(
      colnames(X_cov_net)
    )
    
    return(
      list(
        source = "saved_linear_SCAR_competing_geometry",
        X_cov = X_cov,
        Z_net = Z_net,
        X_cov_net = X_cov_net,
        A_network = A_network,
        network_summary = geometry$network_summary,
        modularity_q = modularity_q,
        modularity_selected_eigenvalues =
          modularity_selected_eigenvalues,
        seeds = seed_row,
        covariate_profile =
          geometry$covariate_profile
      )
    )
  }
  
  make_generated_feature_bundle(
    seed_row = seed_row,
    source_label =
      "regenerated_linear_SCAR_geometry_from_competing_seeds"
  )
}

make_model_specs <- function(bundle) {
  X_cov <- as.matrix(bundle$X_cov)
  X_cov_net <- as.matrix(bundle$X_cov_net)
  
  p_cov <- ncol(X_cov)
  q_net <- ncol(bundle$Z_net)
  
  if (p_cov != 10L) {
    stop(
      "Expected exactly 10 generated covariates; found ",
      p_cov,
      "."
    )
  }
  
  if (q_net < 1L) {
    stop(
      "The network feature matrix must contain at least one column."
    )
  }
  
  groups_cov <- list(
    cov = seq_len(p_cov)
  )
  
  groups_net <- list(
    cov = seq_len(p_cov),
    network = p_cov + seq_len(q_net)
  )
  
  validate_groups(
    groups_cov,
    ncol(X_cov)
  )
  
  validate_groups(
    groups_net,
    ncol(X_cov_net)
  )
  
  geometry_scales <- geometry_half_t_scales(X_cov, bundle$Z_net)

  list(
    cov = list(
      key = "cov",
      method = "Bayesian Linear Cov",
      W = X_cov,
      groups = groups_cov,
      tau_scale = c(cov = TAU_SCALE_COV_SINGLE_ACTIVE),
      tau2_default = c(cov = TAU2_COV_SINGLE_DEFAULT),
      use_pu = FALSE
    ),
    
    cov_pu = list(
      key = "cov_pu",
      method = "Bayesian Linear Cov + PU",
      W = X_cov,
      groups = groups_cov,
      tau_scale = c(cov = TAU_SCALE_COV_SINGLE_ACTIVE),
      tau2_default = c(cov = TAU2_COV_SINGLE_DEFAULT),
      use_pu = TRUE
    ),
    
    cov_net = list(
      key = "cov_net",
      method = "Bayesian Linear Cov + Net",
      W = X_cov_net,
      groups = groups_net,
      tau_scale =
        geometry_scales$tau_scale[c("cov", "network")],
      tau2_default =
        geometry_scales$tau2_default[c("cov", "network")],
      use_pu = FALSE
    ),
    
    cov_net_pu = list(
      key = "cov_net_pu",
      method = "Bayesian Linear Cov + Net + PU",
      W = X_cov_net,
      groups = groups_net,
      tau_scale =
        geometry_scales$tau_scale[c("cov", "network")],
      tau2_default =
        geometry_scales$tau2_default[c("cov", "network")],
      use_pu = TRUE
    )
  )
}

MODEL_REGISTRY <- data.frame(
  model_key = c(
    "cov",
    "cov_pu",
    "cov_net",
    "cov_net_pu"
  ),
  
  method = c(
    "Bayesian Linear Cov",
    "Bayesian Linear Cov + PU",
    "Bayesian Linear Cov + Net",
    "Bayesian Linear Cov + Net + PU"
  ),
  
  stringsAsFactors = FALSE
)

## ============================================================================
## 8. POLYA-GAMMA AND PU HELPERS
## ============================================================================

draw_pg1 <- function(psi) {
  psi <- as.numeric(psi)
  n <- length(psi)
  out <- if (identical(PG_BACKEND, "BayesLogit")) {
    BayesLogit::rpg(num = n, h = rep(1, n), z = psi)
  } else {
    pgdraw::pgdraw(b = 1L, c = psi)
  }
  out <- as.numeric(out)
  if (length(out) != n || any(!is.finite(out)) || any(out <= 0)) {
    stop("The Polya-Gamma backend returned invalid draws.")
  }
  pmax(out, PG_OMEGA_FLOOR)
}

sample_latent_T <- function(Y, psi, eta) {
  Y <- as.numeric(Y)
  p <- plogis(as.numeric(psi))
  eta <- pmin(pmax(as.numeric(eta), ETA_EPS), 1 - ETA_EPS)
  T <- as.integer(Y)
  zero_idx <- which(Y == 0)
  if (length(zero_idx) > 0L) {
    denom <- 1 - (1 - eta) * p[zero_idx]
    denom <- pmax(denom, .Machine$double.eps)
    q_hidden <- eta * p[zero_idx] / denom
    q_hidden <- pmin(pmax(q_hidden, 0), 1)
    T[zero_idx] <- rbinom(length(zero_idx), 1L, q_hidden)
  }
  T[Y == 1] <- 1L
  T
}

sample_eta_given_T <- function(T, Y, a_eta, b_eta) {
  zero_idx <- which(Y == 0)
  n_hidden <- if (length(zero_idx) > 0L) sum(T[zero_idx]) else 0L
  n_observed_positive <- sum(Y == 1)
  eta <- rbeta(
    1L,
    shape1 = a_eta + n_hidden,
    shape2 = b_eta + n_observed_positive
  )
  pmin(pmax(eta, ETA_EPS), 1 - ETA_EPS)
}

## ============================================================================
## 9. GENERIC POLYA-GAMMA GIBBS CHAIN
## ============================================================================

run_pg_chain <- function(
    Y,
    W,
    groups,
    tau_scale,
    tau2_default,
    m_beta,
    s_beta,
    use_pu,
    a_eta,
    b_eta,
    n_burnin,
    n_sample = 0L,
    thin = 1L,
    init_state = NULL,
    initial_eta = INITIAL_ETA,
    nu_tau = NU_TAU,
    state_tail_length = 0L,
    predictive_seed = NULL,
    phase_label = "Bayesian linear",
    burnin_progress_every = 5000L,
    sample_progress_every = 5000L,
    verbose = TRUE
) {
  Y <- as.numeric(Y)
  W <- as.matrix(W)
  n <- length(Y)
  p_features <- ncol(W)
  if (!all(Y %in% c(0, 1))) stop("Y must contain only 0 and 1.")
  if (nrow(W) != n) stop("W and Y have incompatible dimensions.")
  if (p_features < 1L) stop("W must contain at least one feature.")
  validate_groups(groups, p_features)
  group_names <- names(groups)
  tau_scale <- normalize_named_positive(tau_scale, group_names, "tau_scale")
  tau2_default <- normalize_named_positive(
    tau2_default,
    group_names,
    "tau2_default"
  )
  n_burnin <- as.integer(n_burnin)
  n_sample <- as.integer(n_sample)
  thin <- as.integer(thin)
  state_tail_length <- as.integer(state_tail_length)
  if (n_burnin < 0L || n_sample < 0L || thin <= 0L) {
    stop("Invalid MCMC lengths.")
  }
  if (n_sample > 0L && n_sample %% thin != 0L) {
    stop("n_sample must be divisible by thin.")
  }
  total_iterations <- n_burnin + n_sample
  if (state_tail_length < 0L || state_tail_length > total_iterations) {
    stop("state_tail_length is incompatible with the total chain length.")
  }
  Z <- cbind(beta0 = 1, W)
  feature_names <- colnames(W)
  if (is.null(feature_names)) {
    feature_names <- paste0("feature", seq_len(p_features))
  }
  colnames(Z) <- c("beta0", paste0("coef_", feature_names))
  p_beta <- ncol(Z)
  prior_mean <- c(m_beta, rep(0, p_features))
  names(prior_mean) <- colnames(Z)
  if (is.null(init_state)) {
    beta <- prior_mean
    tau2 <- tau2_default
    aux_tau <- setNames(numeric(length(group_names)), group_names)
    for (gname in group_names) {
      aux_tau[gname] <- draw_inverse_gamma_one(
        shape = 1 / 2,
        scale = 1 / (tau_scale[gname]^2)
      )
    }
    if (isTRUE(use_pu)) {
      eta <- pmin(pmax(initial_eta, ETA_EPS), 1 - ETA_EPS)
      T <- as.integer(Y)
    } else {
      eta <- NA_real_
      T <- NULL
    }
  } else {
    beta <- as.numeric(init_state$beta)
    if (length(beta) != p_beta || any(!is.finite(beta))) {
      stop("init_state$beta is incompatible with the current design.")
    }
    names(beta) <- colnames(Z)
    tau2 <- normalize_named_positive(
      init_state$tau2,
      group_names,
      "init_state$tau2"
    )
    aux_tau <- normalize_named_positive(
      init_state$aux_tau,
      group_names,
      "init_state$aux_tau"
    )
    if (isTRUE(use_pu)) {
      eta <- as.numeric(init_state$eta)
      if (length(eta) != 1L || !is.finite(eta)) {
        stop("init_state$eta must be one finite value.")
      }
      eta <- pmin(pmax(eta, ETA_EPS), 1 - ETA_EPS)
      T <- as.integer(init_state$T)
      if (length(T) != n || any(!T %in% c(0L, 1L))) {
        stop("init_state$T must be a binary vector of length n.")
      }
      T[Y == 1] <- 1L
    } else {
      eta <- NA_real_
      T <- NULL
    }
  }
  make_prior_precision <- function(tau2_current) {
    prior_precision <- numeric(p_beta)
    prior_precision[1L] <- 1 / (s_beta^2)
    for (gname in group_names) {
      beta_idx <- 1L + groups[[gname]]
      prior_precision[beta_idx] <- 1 / tau2_current[gname]
    }
    if (any(!is.finite(prior_precision)) || any(prior_precision <= 0)) {
      stop("Invalid Gaussian prior precision.")
    }
    prior_precision
  }
  update_tau_hierarchy <- function(beta_current, tau2_current, aux_current) {
    beta_features <- beta_current[-1L]
    for (gname in group_names) {
      idx <- groups[[gname]]
      tau2_current[gname] <- draw_inverse_gamma_one(
        shape = (nu_tau + length(idx)) / 2,
        scale = (nu_tau / aux_current[gname]) +
          sum(beta_features[idx]^2) / 2
      )
      aux_current[gname] <- draw_inverse_gamma_one(
        shape = (nu_tau + 1) / 2,
        scale = (1 / tau_scale[gname]^2) +
          (nu_tau / tau2_current[gname])
      )
    }
    list(tau2 = tau2_current, aux_tau = aux_current)
  }
  n_saved_expected <- if (n_sample > 0L) n_sample %/% thin else 0L
  score_sum <- numeric(n)
  n_saved <- 0L
  
  eta_draws <- if (isTRUE(use_pu) && n_saved_expected > 0L) {
    numeric(n_saved_expected)
  } else {
    numeric(0)
  }
  
  max_probability_draws <- if (n_saved_expected > 0L) {
    numeric(n_saved_expected)
  } else {
    numeric(0)
  }
  
  q99_probability_draws <- max_probability_draws
  top5_mean_probability_draws <- max_probability_draws
  beta0_draws <- max_probability_draws
  
  tau_draws <- if (n_saved_expected > 0L) {
    matrix(
      NA_real_,
      nrow = n_saved_expected,
      ncol = length(group_names),
      dimnames = list(NULL, paste0("tau_", group_names))
    )
  } else {
    matrix(numeric(0), nrow = 0L, ncol = length(group_names))
  }
  
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
        "predictive_seed is required for posterior predictive labels."
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
  
  tail_beta <- if (state_tail_length > 0L) {
    matrix(
      NA_real_,
      nrow = state_tail_length,
      ncol = p_beta,
      dimnames = list(NULL, colnames(Z))
    )
  } else {
    NULL
  }
  tail_tau2 <- if (state_tail_length > 0L) {
    matrix(
      NA_real_,
      nrow = state_tail_length,
      ncol = length(group_names),
      dimnames = list(NULL, group_names)
    )
  } else {
    NULL
  }
  tail_aux_tau <- if (state_tail_length > 0L) {
    matrix(
      NA_real_,
      nrow = state_tail_length,
      ncol = length(group_names),
      dimnames = list(NULL, group_names)
    )
  } else {
    NULL
  }
  tail_eta <- if (state_tail_length > 0L && isTRUE(use_pu)) {
    numeric(state_tail_length)
  } else {
    numeric(0)
  }
  tail_start_iteration <- if (state_tail_length > 0L) {
    total_iterations - state_tail_length + 1L
  } else {
    Inf
  }
  tail_index <- 0L
  start_time <- Sys.time()
  if (verbose) {
    cat(sprintf("\n[%s] %s sampler started.\n", timestamp_now(), phase_label))
    cat("PU model            :", isTRUE(use_pu), "\n")
    cat("Observations        :", n, "\n")
    cat("Regression columns  :", p_beta, "\n")
    cat("Burn-in iterations  :", n_burnin, "\n")
    cat("Sampling iterations :", n_sample, "\n")
    cat("Thinning            :", thin, "\n")
    cat("Retained draws      :", n_saved_expected, "\n")
    cat("Saved state tail    :", state_tail_length, "\n")
  }
  if (total_iterations > 0L) {
    for (iter in seq_len(total_iterations)) {
      response <- if (isTRUE(use_pu)) T else Y
      psi <- as.vector(Z %*% beta)
      omega <- draw_pg1(psi)
      kappa <- response - 0.5
      prior_precision <- make_prior_precision(tau2)
      weighted_Z <- Z * omega
      Q <- symmetrize(crossprod(Z, weighted_Z)) +
        diag(prior_precision, nrow = p_beta, ncol = p_beta)
      h <- as.vector(crossprod(Z, kappa)) +
        prior_precision * prior_mean
      if (any(!is.finite(Q)) || any(!is.finite(h))) {
        stop("Non-finite Gaussian full conditional for beta.")
      }
      R <- chol_safe(Q)$R
      beta_mean <- as.vector(backsolve(R, forwardsolve(t(R), h)) )
      beta <- beta_mean + as.vector(backsolve(R, rnorm(p_beta)))
      names(beta) <- colnames(Z)
      tau_update <- update_tau_hierarchy(beta, tau2, aux_tau)
      tau2 <- tau_update$tau2
      aux_tau <- tau_update$aux_tau
      if (isTRUE(use_pu)) {
        eta <- sample_eta_given_T(
          T = T,
          Y = Y,
          a_eta = a_eta,
          b_eta = b_eta
        )
        psi_new <- as.vector(Z %*% beta)
        T <- sample_latent_T(Y = Y, psi = psi_new, eta = eta)
      }
      if (state_tail_length > 0L && iter >= tail_start_iteration) {
        tail_index <- tail_index + 1L
        tail_beta[tail_index, ] <- beta
        tail_tau2[tail_index, ] <- tau2[group_names]
        tail_aux_tau[tail_index, ] <- aux_tau[group_names]
        if (isTRUE(use_pu)) tail_eta[tail_index] <- eta
      }
      if (iter > n_burnin) {
        sample_iter <- iter - n_burnin
        if (sample_iter %% thin == 0L) {
          n_saved <- n_saved + 1L
          p_current <- plogis(as.vector(Z %*% beta))
          score_sum <- score_sum + p_current
          max_probability_draws[n_saved] <- max(p_current)
          q99_probability_draws[n_saved] <- unname(
            quantile(p_current, 0.99, names = FALSE)
          )
          top5_mean_probability_draws[n_saved] <- mean(
            sort(p_current, decreasing = TRUE)[1:5]
          )
          beta0_draws[n_saved] <- beta[1L]
          tau_draws[n_saved, ] <- sqrt(tau2[group_names])
          ## Retain eta only as a diagnostic for PU models.
          if (isTRUE(use_pu)) {
            eta_draws[n_saved] <- eta
          }
          
          ## For every model and every unit, construct predictive labels from
          ## T_i ~ Bernoulli(plogis(f_i)).
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
        }
        
        if (verbose && sample_progress_every > 0L &&
            sample_iter %% sample_progress_every == 0L) {
          progress_message(
            paste0(phase_label, " posterior sampling"),
            sample_iter,
            n_sample
          )
        }
      } else if (verbose && burnin_progress_every > 0L &&
                 iter %% burnin_progress_every == 0L) {
        progress_message(
          paste0(phase_label, " burn-in"),
          iter,
          n_burnin
        )
      }
    }
  }
  if (n_saved != n_saved_expected) {
    stop("Retained-draw count mismatch.")
  }
  if (state_tail_length > 0L && tail_index != state_tail_length) {
    stop("Pilot tail was not filled exactly.")
  }
  end_time <- Sys.time()
  runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
  if (verbose) {
    cat(sprintf("[%s] %s sampler finished.\n", timestamp_now(), phase_label))
    cat("Runtime:", format_duration(runtime_sec), "\n")
    cat("Final tau values:\n")
    print(sqrt(tau2))
    if (isTRUE(use_pu)) {
      cat("Final eta:", eta, "\n")
      cat("Final hidden T count:", sum(T[Y == 0]), "\n")
    }
  }
  list(
    state = list(
      beta = beta,
      tau2 = tau2,
      aux_tau = aux_tau,
      eta = if (isTRUE(use_pu)) eta else NULL,
      T = if (isTRUE(use_pu)) T else NULL
    ),
    tail = if (state_tail_length > 0L) {
      list(
        beta = tail_beta,
        tau2 = tail_tau2,
        aux_tau = tail_aux_tau,
        eta = if (isTRUE(use_pu)) tail_eta else NULL
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
    tau_draws = tau_draws,
    n_saved = n_saved,
    time_start = start_time,
    time_end = end_time,
    time_fit_sec = runtime_sec
  )
}

## ============================================================================
## 10. FIVE INDEPENDENT PILOT DATA SETS
## ============================================================================

PILOT_SETTINGS_SIGNATURE <- list(
  script_version = PILOT_CALIBRATION_VERSION,
  master_seed = MASTER_SEED,
  n_pilot = N_PILOT,
  pilot_burnin = PILOT_BURNIN,
  pilot_tail_length = PILOT_TAIL_LENGTH,
  
  data_design = list(
    hidden_covariate_only_idx =
      HIDDEN_COV_ONLY_IDX,
    hidden_network_only_idx =
      HIDDEN_NET_ONLY_IDX,
    hidden_both_idx =
      HIDDEN_BOTH_IDX,
    observed_covariate_only_idx =
      OBSERVED_COV_ONLY_IDX,
    observed_network_only_idx =
      OBSERVED_NET_ONLY_IDX,
    observed_both_idx =
      OBSERVED_BOTH_IDX,
    covariate_signal_idx =
      COVARIATE_SIGNAL_IDX,
    network_signal_idx =
      NETWORK_SIGNAL_IDX,
    observed_y = OBSERVED_Y,
    true_t = TRUE_T
  ),
  
  dgp = list(
    covariates = list(
      type = "ten independent Gaussian covariates",
      n_covariates = N_COVARIATES,
      background_mean =
        COVARIATE_BACKGROUND_MEAN,
      signal_mean =
        COVARIATE_SIGNAL_MEAN,
      common_sd = COMMON_SD,
      allocation =
        "covariate-only plus both-signal positives receive the mean shift"
    ),
    
    network = list(
      type = "two-block Bernoulli SBM",
      block_zero_size =
        N_NETWORK_BACKGROUND,
      block_one_size =
        N_NETWORK_SIGNAL,
      block_definition =
        "network-only plus both-signal positives form block 1",
      p_00 = P_00,
      p_01 = P_01,
      p_11 = P_11,
      triad_closure_prob =
        TRIAD_CLOSURE_PROB,
      triad_max_new_edges =
        TRIAD_MAX_NEW_EDGES
    ),
    
    modularity = list(
      eig_tol = MOD_EIG_TOL,
      max_q = MOD_MAX_Q,
      k_max = MOD_K_MAX
    )
  ),
  
  priors = list(
    beta0 = beta0_prior_for_sampler,
    eta = eta_prior_for_sampler,
    tau_scale_rule = TAU_SCALE_RULE,
    nu_tau = NU_TAU
  ),
  
  pg_backend = PG_BACKEND
)

run_or_load_one_pilot_model <- function(pilot_number, bundle, spec) {
  checkpoint_file <- file.path(
    PILOT_CHECKPOINT_DIR,
    sprintf("pilot_%02d_%s.rds", pilot_number, spec$key)
  )
  if (isTRUE(RESUME_PILOTS) && !isTRUE(FORCE_FRESH_PILOTS) &&
      file.exists(checkpoint_file)) {
    candidate <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
    if (is.list(candidate) &&
        identical(candidate$settings_signature, PILOT_SETTINGS_SIGNATURE) &&
        identical(as.integer(candidate$pilot), as.integer(pilot_number)) &&
        identical(candidate$model_key, spec$key)) {
      cat(sprintf(
        "[%s] Reusing pilot %02d | %s\n",
        timestamp_now(), pilot_number, spec$method
      ))
      return(candidate)
    }
  }
  mcmc_seed <- as.integer(
    PILOT_MCMC_SEED_BASE +
      MODEL_SEED_OFFSETS[spec$key] +
      10000L * pilot_number
  )
  set.seed(mcmc_seed)
  fit <- run_pg_chain(
    Y = OBSERVED_Y,
    W = spec$W,
    groups = spec$groups,
    tau_scale = spec$tau_scale,
    tau2_default = spec$tau2_default,
    m_beta = beta0_prior_for_sampler$mean,
    s_beta = beta0_prior_for_sampler$sd,
    use_pu = spec$use_pu,
    a_eta = eta_prior_for_sampler$a_eta,
    b_eta = eta_prior_for_sampler$b_eta,
    n_burnin = PILOT_BURNIN,
    n_sample = 0L,
    thin = 1L,
    init_state = NULL,
    initial_eta = INITIAL_ETA,
    nu_tau = NU_TAU,
    state_tail_length = PILOT_TAIL_LENGTH,
    phase_label = paste0(
      "Pilot ", sprintf("%02d", pilot_number), " | ", spec$method
    ),
    burnin_progress_every = PILOT_PROGRESS_EVERY,
    sample_progress_every = 0L,
    verbose = TRUE
  )
  result <- list(
    settings_signature = PILOT_SETTINGS_SIGNATURE,
    pilot = pilot_number,
    model_key = spec$key,
    method = spec$method,
    feature_count = ncol(spec$W),
    network_feature_count = ncol(bundle$Z_net),
    mcmc_seed = mcmc_seed,
    tail = fit$tail,
    final_state = fit$state,
    runtime_sec = fit$time_fit_sec
  )
  saveRDS(result, checkpoint_file)
  result
}

run_all_pilots <- function() {
  pilot_results <- vector("list", N_PILOT)
  for (pp in seq_len(N_PILOT)) {
    seed_row <- pilot_seed_table[
      pilot_seed_table$pilot == pp,
      ,
      drop = FALSE
    ]
    cat("\n\n####################################################################\n")
    cat(sprintf("# GLOBAL PILOT DATA SET %d / %d\n", pp, N_PILOT))
    cat("####################################################################\n")
    bundle <- make_generated_feature_bundle(
      seed_row = seed_row,
      source_label = "independent_global_pilot"
    )
    if (isTRUE(SAVE_PILOT_GEOMETRIES)) {
      saveRDS(
        list(
          pilot = pp,
          seeds = seed_row,
          T = TRUE_T,
          Y = OBSERVED_Y,
          hidden_positive_idx = HIDDEN_POS_IDX,
          observed_positive_idx = OBSERVED_POS_IDX,
          covariate_signal_idx = COVARIATE_SIGNAL_IDX,
          network_signal_idx = NETWORK_SIGNAL_IDX,
          network_block = NETWORK_BLOCK,
          covariate_dgp = PILOT_SETTINGS_SIGNATURE$dgp$covariates,
          network_dgp = PILOT_SETTINGS_SIGNATURE$dgp$network,
          X_cov = bundle$X_cov,
          Z_net = bundle$Z_net,
          A_network = bundle$A_network,
          network_summary = bundle$network_summary,
          modularity_q = bundle$modularity_q
        ),
        file.path(PILOT_GEOMETRY_DIR, sprintf("pilot_%02d_geometry.rds", pp))
      )
    }
    specs <- make_model_specs(bundle)
    model_results <- lapply(
      MODEL_REGISTRY$model_key,
      function(key) run_or_load_one_pilot_model(pp, bundle, specs[[key]])
    )
    names(model_results) <- MODEL_REGISTRY$model_key
    pilot_results[[pp]] <- model_results
    rm(bundle, specs, model_results)
    invisible(gc())
  }
  pilot_results
}

build_global_calibration <- function(pilot_results) {
  calibration_by_model <- list()
  summary_rows <- list()
  summary_index <- 1L
  per_pilot_rows <- list()
  per_pilot_index <- 1L
  for (key in MODEL_REGISTRY$model_key) {
    model_tail_objects <- lapply(
      pilot_results,
      function(pilot_result) pilot_result[[key]]$tail
    )
    beta_fixed_pooled <- do.call(
      rbind,
      lapply(model_tail_objects, function(tail_obj) {
        tail_obj$beta[, seq_len(11L), drop = FALSE]
      })
    )
    if (ncol(beta_fixed_pooled) != 11L) {
      stop("Pilot beta pooling failed for model ", key, ".")
    }
    group_names <- colnames(model_tail_objects[[1L]]$tau2)
    tau2_pooled <- do.call(
      rbind,
      lapply(model_tail_objects, function(tail_obj) {
        tail_obj$tau2[, group_names, drop = FALSE]
      })
    )
    aux_pooled <- do.call(
      rbind,
      lapply(model_tail_objects, function(tail_obj) {
        tail_obj$aux_tau[, group_names, drop = FALSE]
      })
    )
    eta_pooled <- if (!is.null(model_tail_objects[[1L]]$eta)) {
      unlist(lapply(model_tail_objects, `[[`, "eta"), use.names = FALSE)
    } else {
      NULL
    }
    global_now <- list(
      beta0 = median(beta_fixed_pooled[, 1L]),
      cov_beta = apply(beta_fixed_pooled[, 2:11, drop = FALSE], 2L, median),
      tau2 = setNames(apply(tau2_pooled, 2L, median), group_names),
      aux_tau = setNames(apply(aux_pooled, 2L, median), group_names),
      eta = if (!is.null(eta_pooled)) median(eta_pooled) else NULL,
      n_pilot = N_PILOT,
      pooled_tail_draws = nrow(beta_fixed_pooled)
    )
    calibration_by_model[[key]] <- global_now
    parameter_vectors <- list(
      beta0 = beta_fixed_pooled[, 1L],
      tau2_cov = tau2_pooled[, "cov"],
      aux_tau_cov = aux_pooled[, "cov"]
    )
    for (jj in seq_len(10L)) {
      parameter_vectors[[paste0("cov_beta_", jj)]] <- beta_fixed_pooled[, jj + 1L]
    }
    if ("network" %in% group_names) {
      parameter_vectors$tau2_network <- tau2_pooled[, "network"]
      parameter_vectors$aux_tau_network <- aux_pooled[, "network"]
    }
    if (!is.null(eta_pooled)) parameter_vectors$eta <- eta_pooled
    for (parameter_name in names(parameter_vectors)) {
      values <- parameter_vectors[[parameter_name]]
      stats_now <- summary_stats(values)
      summary_rows[[summary_index]] <- data.frame(
        model_key = key,
        method = MODEL_REGISTRY$method[
          match(key, MODEL_REGISTRY$model_key)
        ],
        parameter = parameter_name,
        mean = unname(stats_now["mean"]),
        sd = unname(stats_now["sd"]),
        q025 = unname(stats_now["q025"]),
        median = unname(stats_now["median"]),
        q975 = unname(stats_now["q975"]),
        pooled_tail_draws = length(values),
        stringsAsFactors = FALSE
      )
      summary_index <- summary_index + 1L
    }
    for (pp in seq_len(N_PILOT)) {
      tail_obj <- model_tail_objects[[pp]]
      row_now <- data.frame(
        pilot = pp,
        model_key = key,
        method = MODEL_REGISTRY$method[
          match(key, MODEL_REGISTRY$model_key)
        ],
        beta0_median = median(tail_obj$beta[, 1L]),
        eta_median = if (!is.null(tail_obj$eta)) median(tail_obj$eta) else NA_real_,
        tau2_cov_median = median(tail_obj$tau2[, "cov"]),
        aux_tau_cov_median = median(tail_obj$aux_tau[, "cov"]),
        tau2_network_median = if ("network" %in% colnames(tail_obj$tau2)) {
          median(tail_obj$tau2[, "network"])
        } else {
          NA_real_
        },
        aux_tau_network_median = if ("network" %in% colnames(tail_obj$aux_tau)) {
          median(tail_obj$aux_tau[, "network"])
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
      per_pilot_rows[[per_pilot_index]] <- row_now
      per_pilot_index <- per_pilot_index + 1L
    }
  }
  list(
    settings_signature = PILOT_SETTINGS_SIGNATURE,
    created_at = timestamp_now(),
    calibration_by_model = calibration_by_model,
    distribution_summary = do.call(rbind, summary_rows),
    per_pilot_summary = do.call(rbind, per_pilot_rows),
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
      cat("\nLoaded the completed five-pilot global calibration.\n")
      return(candidate)
    }
  }
  pilot_results <- run_all_pilots()
  calibration <- build_global_calibration(pilot_results)
  saveRDS(calibration, GLOBAL_CALIBRATION_FILE)
  global_pilot_calibration <- calibration
  save(
    global_pilot_calibration,
    file = GLOBAL_CALIBRATION_RDATA
  )
  write.csv(
    calibration$distribution_summary,
    file.path(OUT_DIR, "global_pilot_calibration_distribution_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    calibration$per_pilot_summary,
    file.path(OUT_DIR, "global_pilot_calibration_per_pilot_summary.csv"),
    row.names = FALSE
  )
  rm(pilot_results)
  invisible(gc())
  calibration
}

global_pilot_calibration <- load_or_build_global_calibration()
GLOBAL_CALIBRATION_HASH <- file_md5_or_na(GLOBAL_CALIBRATION_FILE)

## Ensure the human-readable calibration summaries are present even when the
## completed RDS calibration is reused from an earlier execution.
write.csv(
  global_pilot_calibration$distribution_summary,
  file.path(OUT_DIR, "global_pilot_calibration_distribution_summary.csv"),
  row.names = FALSE
)
write.csv(
  global_pilot_calibration$per_pilot_summary,
  file.path(OUT_DIR, "global_pilot_calibration_per_pilot_summary.csv"),
  row.names = FALSE
)

cat("\n================ GLOBAL PILOT INITIAL VALUES ================\n")
for (key in MODEL_REGISTRY$model_key) {
  cat("\n", MODEL_REGISTRY$method[match(key, MODEL_REGISTRY$model_key)], "\n", sep = "")
  print(global_pilot_calibration$calibration_by_model[[key]])
}

## ============================================================================
## 11. PRODUCTION INITIAL STATES
## ============================================================================

make_global_initial_state <- function(spec, calibration, Y) {
  p_cov <- 10L
  p_features <- ncol(spec$W)
  n_network_features <- p_features - p_cov
  if (n_network_features < 0L) {
    stop("Current production design has fewer than ten covariate columns.")
  }
  beta <- c(
    calibration$beta0,
    as.numeric(calibration$cov_beta),
    if (n_network_features > 0L) rep(0, n_network_features) else numeric(0)
  )
  if (length(beta) != 1L + p_features) {
    stop("Global beta initialization has the wrong length for ", spec$method, ".")
  }
  group_names <- names(spec$groups)
  tau2 <- normalize_named_positive(
    calibration$tau2,
    group_names,
    "global calibration tau2"
  )
  aux_tau <- normalize_named_positive(
    calibration$aux_tau,
    group_names,
    "global calibration aux_tau"
  )
  if (isTRUE(spec$use_pu)) {
    eta <- as.numeric(calibration$eta)
    if (length(eta) != 1L || !is.finite(eta)) {
      stop("Global eta initialization is invalid for ", spec$method, ".")
    }
    ## No latent state is transferred from any pilot.
    T <- as.integer(Y)
  } else {
    eta <- NULL
    T <- NULL
  }
  list(
    beta = beta,
    tau2 = tau2,
    aux_tau = aux_tau,
    eta = eta,
    T = T
  )
}

## ============================================================================
## 12. COMMON EVALUATION
## ============================================================================

evaluate_bayesian_method <- function(
    fit,
    spec,
    replicate_number,
    runtime_sec
) {
  probabilities <- validate_probability_vector(
    fit$posterior_mean_scores,
    N_NODES,
    spec$method
  )
  
  predictive_T_draws <- validate_predictive_label_draws(
    predictive_T_draws = fit$predictive_T_draws,
    n_draws_expected = fit$n_saved,
    n_nodes_expected = N_NODES,
    method_label = spec$method
  )
  
  unlabeled_idx <- which(OBSERVED_Y == 0L)
  hidden_truth_unlabeled <- as.integer(
    TRUE_T[unlabeled_idx] == 1L
  )
  
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
    hidden_idx = HIDDEN_POS_IDX,
    level = 0.80
  )
  
  predictive_95 <- predictive_label_metrics(
    predictive_T_draws = predictive_T_draws,
    truth = TRUE_T,
    hidden_idx = HIDDEN_POS_IDX,
    level = 0.95
  )
  
  eta_mean <- if (
    isTRUE(spec$use_pu) &&
    length(fit$eta_draws) > 0L
  ) {
    mean(fit$eta_draws)
  } else {
    NA_real_
  }
  
  metrics <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    model_key = spec$key,
    method = spec$method,
    
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
      HIDDEN_POS_IDX
    ),
    
    Overall_Brier = overall_brier(
      probabilities,
      TRUE_T
    ),
    
    Hidden_Brier = hidden_brier(
      probabilities,
      HIDDEN_POS_IDX
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
    
    eta_posterior_mean = eta_mean,
    
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
    runtime_sec = runtime_sec,
    stringsAsFactors = FALSE
  )
  
  node_scores <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    model_key = spec$key,
    method = spec$method,
    node = seq_len(N_NODES),
    T = TRUE_T,
    Y = OBSERVED_Y,
    is_hidden_positive = as.integer(
      seq_len(N_NODES) %in% HIDDEN_POS_IDX
    ),
    covariate_signal = as.integer(
      seq_len(N_NODES) %in% COVARIATE_SIGNAL_IDX
    ),
    network_signal = as.integer(
      seq_len(N_NODES) %in% NETWORK_SIGNAL_IDX
    ),
    signal_type = SIGNAL_TYPE,
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
## 13. PRODUCTION CHECKPOINT SIGNATURE
## ============================================================================

PRODUCTION_SETTINGS_SIGNATURE <- list(
  script_version = SCRIPT_VERSION,
  master_seed = MASTER_SEED,
  n_rep = N_REP,
  production_burnin = PRODUCTION_BURNIN,
  production_sampling = PRODUCTION_SAMPLING,
  production_thin = PRODUCTION_THIN,
  global_calibration_hash = GLOBAL_CALIBRATION_HASH,
  use_saved_competing_geometries = USE_SAVED_COMPETING_GEOMETRIES,
  data_design = list(
    true_t = TRUE_T,
    observed_y = OBSERVED_Y,
    hidden_covariate_only_idx = HIDDEN_COV_ONLY_IDX,
    hidden_network_only_idx = HIDDEN_NET_ONLY_IDX,
    hidden_both_idx = HIDDEN_BOTH_IDX,
    observed_covariate_only_idx = OBSERVED_COV_ONLY_IDX,
    observed_network_only_idx = OBSERVED_NET_ONLY_IDX,
    observed_both_idx = OBSERVED_BOTH_IDX,
    covariate_signal_idx = COVARIATE_SIGNAL_IDX,
    network_signal_idx = NETWORK_SIGNAL_IDX,
    network_block = NETWORK_BLOCK
  ),
  dgp = PILOT_SETTINGS_SIGNATURE$dgp,
  priors = list(
    beta0 = beta0_prior_for_sampler,
    eta = eta_prior_for_sampler,
    tau_scale_rule = TAU_SCALE_RULE,
    nu_tau = NU_TAU
  ),
  pg_backend = PG_BACKEND,
  evaluation = list(
    hidden_top_k = HIDDEN_TOP_K,
    overall_recall_k = OVERALL_RECALL_K,
    predictive_interval_levels = PREDICTIVE_INTERVAL_LEVELS,
    predictive_quantile_type = PREDICTIVE_QUANTILE_TYPE,
    truth_for_coverage =
      "known latent TRUE_T; hidden positives count as T=1",
    predictive_label_draws = paste0(
      "For every model, unit and retained posterior beta draw, ",
      "sample T_i from Bernoulli(plogis(f_i)), with f_i = Z_i beta."
    ),
    H_interval = "not computed",
    poisson_binomial = "not used"
  )
)

## ============================================================================
## 14. ONE COMPLETE PRODUCTION REPLICATION
## ============================================================================

run_one_production_replication <- function(replicate_number) {
  checkpoint_file <- file.path(
    PRODUCTION_CHECKPOINT_DIR,
    sprintf("replicate_%03d.rds", replicate_number)
  )
  checkpoint <- NULL
  if (isTRUE(RESUME_PRODUCTION) && !isTRUE(FORCE_FRESH_PRODUCTION) &&
      file.exists(checkpoint_file)) {
    candidate <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
    if (is.list(candidate) &&
        identical(candidate$settings_signature, PRODUCTION_SETTINGS_SIGNATURE) &&
        identical(as.integer(candidate$replicate), as.integer(replicate_number))) {
      checkpoint <- candidate
    }
  }
  if (!is.null(checkpoint) && isTRUE(checkpoint$complete)) {
    return(checkpoint)
  }
  replication_start <- Sys.time()
  bundle <- load_or_generate_production_bundle(replicate_number)
  specs <- make_model_specs(bundle)
  if (is.null(checkpoint)) {
    checkpoint <- list(
      settings_signature = PRODUCTION_SETTINGS_SIGNATURE,
      replicate = replicate_number,
      complete = FALSE,
      method_results = list(),
      method_failures = list(),
      geometry_summary = NULL,
      modularity_summary = NULL,
      time_start = replication_start,
      time_end = NULL,
      runtime_sec = NA_real_
    )
  }
  network_summary <- bundle$network_summary
  if (is.null(network_summary)) {
    network_summary <- data.frame(stringsAsFactors = FALSE)
  } else {
    network_summary <- as.data.frame(network_summary, stringsAsFactors = FALSE)
  }
  geometry_base <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    geometry_source = bundle$source,
    covariate_seed = bundle$seeds$covariate_seed,
    network_seed = bundle$seeds$network_seed,
    modularity_seed = bundle$seeds$modularity_seed,
    n = N_NODES,
    n_true_zero = N_TRUE_ZERO,
    n_true_positive = N_TRUE_POS,
    n_observed_positive = N_OBSERVED_POS,
    n_hidden_positive = N_HIDDEN,
    n_covariate_signal = length(COVARIATE_SIGNAL_IDX),
    n_network_signal = length(NETWORK_SIGNAL_IDX),
    n_both_signal = length(intersect(COVARIATE_SIGNAL_IDX, NETWORK_SIGNAL_IDX)),
    stringsAsFactors = FALSE
  )
  checkpoint$geometry_summary <- if (ncol(network_summary) > 0L) {
    cbind(geometry_base, network_summary)
  } else {
    geometry_base
  }
  checkpoint$modularity_summary <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    q_network = bundle$modularity_q,
    selected_eigenvalues = if (
      length(bundle$modularity_selected_eigenvalues) > 0L
    ) {
      paste(signif(bundle$modularity_selected_eigenvalues, 8), collapse = ";")
    } else {
      ""
    },
    n_covariates = ncol(bundle$X_cov),
    n_network_features = ncol(bundle$Z_net),
    n_covariates_plus_network = ncol(bundle$X_cov_net),
    stringsAsFactors = FALSE
  )
  for (key in MODEL_REGISTRY$model_key) {
    spec <- specs[[key]]
    if (key %in% names(checkpoint$method_results)) {
      cat(sprintf(
        "[%s] Replicate %03d/%03d | reusing %s\n",
        timestamp_now(), replicate_number, N_REP, spec$method
      ))
      next
    }
    cat(sprintf(
      "\n[%s] Replicate %03d/%03d | fitting %s\n",
      timestamp_now(), replicate_number, N_REP, spec$method
    ))
    model_start <- Sys.time()
    model_error <- NA_character_
    calibration <- global_pilot_calibration$calibration_by_model[[key]]
    init_state <- make_global_initial_state(
      spec = spec,
      calibration = calibration,
      Y = OBSERVED_Y
    )
    mcmc_seed <- as.integer(
      PRODUCTION_MCMC_SEED_BASE +
        MODEL_SEED_OFFSETS[key] +
        10000L * replicate_number
    )
    predictive_seed <- as.integer(
      PRODUCTION_PREDICTIVE_SEED_BASE +
        MODEL_SEED_OFFSETS[key] +
        10000L * replicate_number
    )
    fit <- tryCatch({
      set.seed(mcmc_seed)
      run_pg_chain(
        Y = OBSERVED_Y,
        W = spec$W,
        groups = spec$groups,
        tau_scale = spec$tau_scale,
        tau2_default = spec$tau2_default,
        m_beta = beta0_prior_for_sampler$mean,
        s_beta = beta0_prior_for_sampler$sd,
        use_pu = spec$use_pu,
        a_eta = eta_prior_for_sampler$a_eta,
        b_eta = eta_prior_for_sampler$b_eta,
        n_burnin = PRODUCTION_BURNIN,
        n_sample = PRODUCTION_SAMPLING,
        thin = PRODUCTION_THIN,
        init_state = init_state,
        initial_eta = INITIAL_ETA,
        nu_tau = NU_TAU,
        state_tail_length = 0L,
        predictive_seed = predictive_seed,
        phase_label = paste0(
          "Replication ", sprintf("%03d", replicate_number),
          " | ", spec$method
        ),
        burnin_progress_every = PRODUCTION_BURNIN_PROGRESS_EVERY,
        sample_progress_every = PRODUCTION_SAMPLE_PROGRESS_EVERY,
        verbose = TRUE
      )
    }, error = function(e) {
      model_error <<- conditionMessage(e)
      NULL
    })
    runtime_sec <- as.numeric(difftime(Sys.time(), model_start, units = "secs"))
    if (is.null(fit)) {
      checkpoint$method_failures[[key]] <- data.frame(
        scenario = SCENARIO_NAME,
        replicate = replicate_number,
        model_key = key,
        method = spec$method,
        runtime_sec = runtime_sec,
        error_message = model_error,
        stringsAsFactors = FALSE
      )
      saveRDS(checkpoint, checkpoint_file)
      if (isTRUE(FAIL_IF_ANY_MODEL_FAILS)) {
        stop(
          "Model failure in replicate ", replicate_number,
          " for ", spec$method, ": ", model_error,
          ". Partial checkpoint saved to ", checkpoint_file, "."
        )
      }
      next
    }
    evaluation <- evaluate_bayesian_method(
      fit = fit,
      spec = spec,
      replicate_number = replicate_number,
      runtime_sec = runtime_sec
    )
    checkpoint$method_results[[key]] <- list(
      metrics = evaluation$metrics,
      node_scores = evaluation$node_scores,
      final_state = fit$state,
      n_saved = fit$n_saved,
      mcmc_seed = mcmc_seed,
      predictive_seed = predictive_seed
    )
    checkpoint$method_failures[[key]] <- NULL
    saveRDS(checkpoint, checkpoint_file)
    rm(fit, evaluation, init_state)
    invisible(gc())
  }
  missing_models <- setdiff(
    MODEL_REGISTRY$model_key,
    names(checkpoint$method_results)
  )
  if (length(missing_models) > 0L) {
    checkpoint$complete <- FALSE
    saveRDS(checkpoint, checkpoint_file)
    stop(
      "Replication ", replicate_number,
      " did not complete all four models. Missing: ",
      paste(missing_models, collapse = ", ")
    )
  }
  checkpoint$complete <- TRUE
  checkpoint$time_end <- Sys.time()
  checkpoint$runtime_sec <- as.numeric(difftime(
    checkpoint$time_end,
    checkpoint$time_start,
    units = "secs"
  ))
  saveRDS(checkpoint, checkpoint_file)
  rm(bundle, specs)
  invisible(gc())
  checkpoint
}

## ============================================================================
## 15. RUN ALL 100 PRODUCTION REPLICATIONS
## ============================================================================

production_results <- vector("list", N_REP)
production_start <- Sys.time()
n_loaded_complete <- 0L
n_computed_now <- 0L

for (rr in seq_len(N_REP)) {
  checkpoint_file <- file.path(
    PRODUCTION_CHECKPOINT_DIR,
    sprintf("replicate_%03d.rds", rr)
  )
  was_complete_before <- FALSE
  if (isTRUE(RESUME_PRODUCTION) && !isTRUE(FORCE_FRESH_PRODUCTION) &&
      file.exists(checkpoint_file)) {
    prior_checkpoint <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
    was_complete_before <- is.list(prior_checkpoint) &&
      identical(prior_checkpoint$settings_signature, PRODUCTION_SETTINGS_SIGNATURE) &&
      isTRUE(prior_checkpoint$complete)
  }
  production_results[[rr]] <- run_one_production_replication(rr)
  if (was_complete_before) {
    n_loaded_complete <- n_loaded_complete + 1L
  } else {
    n_computed_now <- n_computed_now + 1L
  }
  elapsed_sec <- as.numeric(difftime(Sys.time(), production_start, units = "secs"))
  projected_total_sec <- elapsed_sec / rr * N_REP
  cat(sprintf(
    "[%s] Completed production replication %03d/%03d | elapsed %s | projected total %s\n",
    timestamp_now(), rr, N_REP,
    format_duration(elapsed_sec),
    format_duration(projected_total_sec)
  ))
}

production_end <- Sys.time()

## ============================================================================
## 16. COMBINE PRODUCTION RESULTS
## ============================================================================

if (any(!vapply(production_results, function(x) is.list(x) && isTRUE(x$complete), logical(1)))) {
  stop("Not all 100 Bayesian linear production replications completed.")
}

metrics_by_replicate <- do.call(
  rbind,
  lapply(production_results, function(rep_result) {
    do.call(rbind, lapply(MODEL_REGISTRY$model_key, function(key) {
      rep_result$method_results[[key]]$metrics
    }))
  })
)
row.names(metrics_by_replicate) <- NULL

scores_all_nodes <- do.call(
  rbind,
  lapply(production_results, function(rep_result) {
    do.call(rbind, lapply(MODEL_REGISTRY$model_key, function(key) {
      rep_result$method_results[[key]]$node_scores
    }))
  })
)
row.names(scores_all_nodes) <- NULL

geometry_summary <- do.call(rbind, lapply(production_results, `[[`, "geometry_summary"))
modularity_summary <- do.call(rbind, lapply(production_results, `[[`, "modularity_summary"))
row.names(geometry_summary) <- NULL
row.names(modularity_summary) <- NULL

failure_parts <- Filter(
  Negate(is.null),
  lapply(production_results, function(rep_result) {
    if (length(rep_result$method_failures) == 0L) return(NULL)
    do.call(rbind, rep_result$method_failures)
  })
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

expected_metric_rows <- N_REP * nrow(MODEL_REGISTRY)
expected_score_rows <- expected_metric_rows * N_NODES
if (nrow(metrics_by_replicate) != expected_metric_rows) {
  stop("The combined metric table does not contain 100 x 4 rows.")
}
if (nrow(scores_all_nodes) != expected_score_rows) {
  stop("The combined node-score table does not contain 100 x 4 x 400 rows.")
}

method_order <- MODEL_REGISTRY$method
metrics_by_replicate$method <- factor(metrics_by_replicate$method, levels = method_order)
metrics_by_replicate <- metrics_by_replicate[
  order(metrics_by_replicate$method, metrics_by_replicate$replicate),
  ,
  drop = FALSE
]
metrics_by_replicate$method <- as.character(metrics_by_replicate$method)

scores_all_nodes$method <- factor(scores_all_nodes$method, levels = method_order)
scores_all_nodes <- scores_all_nodes[
  order(scores_all_nodes$method, scores_all_nodes$replicate, scores_all_nodes$node),
  ,
  drop = FALSE
]
scores_all_nodes$method <- as.character(scores_all_nodes$method)

write.csv(
  metrics_by_replicate,
  file.path(OUT_DIR, "Bayesian_linear_metrics_by_replicate.csv"),
  row.names = FALSE
)
write.csv(
  scores_all_nodes,
  file.path(OUT_DIR, "Bayesian_linear_all_node_probabilities.csv"),
  row.names = FALSE
)
write.csv(
  geometry_summary,
  file.path(OUT_DIR, "Bayesian_linear_geometry_summary.csv"),
  row.names = FALSE
)
write.csv(
  modularity_summary,
  file.path(OUT_DIR, "Bayesian_linear_modularity_summary.csv"),
  row.names = FALSE
)
write.csv(
  failure_log,
  file.path(OUT_DIR, "Bayesian_linear_failure_log.csv"),
  row.names = FALSE
)

## ============================================================================
## 17. MONTE CARLO SUMMARIES
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

summary_rows <- list()
summary_index <- 1L

for (method_name in method_order) {
  method_data <- metrics_by_replicate[
    metrics_by_replicate$method == method_name,
    ,
    drop = FALSE
  ]
  
  if (nrow(method_data) != N_REP) {
    stop(
      "Method ",
      method_name,
      " does not have exactly 100 replicate rows."
    )
  }
  
  for (metric_name in main_metric_names) {
    statistics <- summary_stats(
      method_data[[metric_name]]
    )
    
    summary_rows[[summary_index]] <- data.frame(
      scenario = SCENARIO_NAME,
      method = method_name,
      metric = metric_name,
      direction = unname(metric_direction[metric_name]),
      n_rep = sum(is.finite(method_data[[metric_name]])),
      mean = unname(statistics["mean"]),
      sd = unname(statistics["sd"]),
      mc_se = unname(statistics["mc_se"]),
      q025 = unname(statistics["q025"]),
      median = unname(statistics["median"]),
      q975 = unname(statistics["q975"]),
      stringsAsFactors = FALSE
    )
    
    summary_index <- summary_index + 1L
  }
}

metric_summary_long <- do.call(
  rbind,
  summary_rows
)
row.names(metric_summary_long) <- NULL

write.csv(
  metric_summary_long,
  file.path(
    OUT_DIR,
    "Bayesian_linear_metric_summary_long.csv"
  ),
  row.names = FALSE
)

main_table_means <- data.frame(
  scenario = SCENARIO_NAME,
  method = method_order,
  stringsAsFactors = FALSE
)

for (metric_name in main_metric_names) {
  metric_rows <- metric_summary_long[
    metric_summary_long$metric == metric_name,
    ,
    drop = FALSE
  ]
  
  metric_rows <- metric_rows[
    match(method_order, metric_rows$method),
    ,
    drop = FALSE
  ]
  
  main_table_means[[metric_name]] <- metric_rows$mean
}

write.csv(
  main_table_means,
  file.path(
    OUT_DIR,
    "Bayesian_linear_FINAL_MAIN_TABLE_means.csv"
  ),
  row.names = FALSE
)

format_mean_sd <- function(
    mean_value,
    sd_value,
    digits = 3L
) {
  if (!is.finite(mean_value)) return("")
  
  if (!is.finite(sd_value)) {
    return(
      formatC(
        mean_value,
        format = "f",
        digits = digits
      )
    )
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
  Overall_Recall_at_45 = "Overall Recall@45",
  Hidden_Recall_at_15 = "Hidden Recall@15",
  Overall_Coverage_80 = "Overall Coverage 80%",
  Hidden_Coverage_80 = "Hidden Coverage 80%",
  Overall_Interval_Score_80 =
    "Overall Interval Score 80%",
  Hidden_Interval_Score_80 =
    "Hidden Interval Score 80%",
  Overall_Coverage_95 = "Overall Coverage 95%",
  Hidden_Coverage_95 = "Hidden Coverage 95%",
  Overall_Interval_Score_95 =
    "Overall Interval Score 95%",
  Hidden_Interval_Score_95 =
    "Hidden Interval Score 95%"
)

paper_table_mean_sd <- data.frame(
  Method = method_order,
  stringsAsFactors = FALSE
)

coverage_metrics <- c(
  "Overall_Coverage_80",
  "Hidden_Coverage_80",
  "Overall_Coverage_95",
  "Hidden_Coverage_95"
)

for (metric_name in main_metric_names) {
  metric_rows <- metric_summary_long[
    metric_summary_long$metric == metric_name,
    ,
    drop = FALSE
  ]
  
  metric_rows <- metric_rows[
    match(method_order, metric_rows$method),
    ,
    drop = FALSE
  ]
  
  if (metric_name %in% coverage_metrics) {
    formatted <- vapply(
      metric_rows$mean,
      function(x) {
        if (!is.finite(x)) {
          ""
        } else {
          formatC(
            x,
            format = "f",
            digits = 3L
          )
        }
      },
      character(1)
    )
  } else {
    formatted <- mapply(
      format_mean_sd,
      metric_rows$mean,
      metric_rows$sd,
      MoreArgs = list(digits = 3L),
      USE.NAMES = FALSE
    )
  }
  
  paper_table_mean_sd[[paper_column_labels[metric_name]]] <- formatted
}

write.csv(
  paper_table_mean_sd,
  file.path(
    OUT_DIR,
    "Bayesian_linear_FINAL_PAPER_TABLE_mean_SD.csv"
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

uncertainty_table_mean_sd <- paper_table_mean_sd[
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
    "Bayesian_linear_uncertainty_metrics_80_95.csv"
  ),
  row.names = FALSE,
  na = ""
)

additional_metric_names <- c(
  "eta_posterior_mean",
  "anchor_max_median",
  "anchor_Pr_over_090",
  "anchor_Pr_over_095",
  "anchor_Pr_over_099",
  "anchor_q99_median",
  "anchor_top5_median",
  "beta0_posterior_mean",
  "runtime_sec"
)

additional_summary_rows <- list()
additional_index <- 1L

for (method_name in method_order) {
  method_data <- metrics_by_replicate[
    metrics_by_replicate$method == method_name,
    ,
    drop = FALSE
  ]
  
  for (metric_name in additional_metric_names) {
    statistics <- summary_stats(
      method_data[[metric_name]]
    )
    
    additional_summary_rows[[additional_index]] <- data.frame(
      method = method_name,
      metric = metric_name,
      n_rep = sum(is.finite(method_data[[metric_name]])),
      mean = unname(statistics["mean"]),
      sd = unname(statistics["sd"]),
      mc_se = unname(statistics["mc_se"]),
      q025 = unname(statistics["q025"]),
      median = unname(statistics["median"]),
      q975 = unname(statistics["q975"]),
      stringsAsFactors = FALSE
    )
    
    additional_index <- additional_index + 1L
  }
}

additional_diagnostic_summary <- do.call(
  rbind,
  additional_summary_rows
)

write.csv(
  additional_diagnostic_summary,
  file.path(
    OUT_DIR,
    "Bayesian_linear_additional_diagnostic_summary.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 18. SAVE COMPLETE OUTPUT BUNDLE
## ============================================================================

total_production_runtime_sec <- as.numeric(difftime(
  production_end,
  production_start,
  units = "secs"
))

run_config <- list(
  script_version = SCRIPT_VERSION,
  pilot_calibration_version = PILOT_CALIBRATION_VERSION,
  created_at = timestamp_now(),
  scenario = SCENARIO_NAME,
  prior_rule = TAU_SCALE_RULE,
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
  observed_labels = OBSERVED_Y,
  signal_type = SIGNAL_TYPE,
  hidden_covariate_only_idx = HIDDEN_COV_ONLY_IDX,
  hidden_network_only_idx = HIDDEN_NET_ONLY_IDX,
  hidden_both_idx = HIDDEN_BOTH_IDX,
  observed_covariate_only_idx = OBSERVED_COV_ONLY_IDX,
  observed_network_only_idx = OBSERVED_NET_ONLY_IDX,
  observed_both_idx = OBSERVED_BOTH_IDX,
  covariate_signal_idx = COVARIATE_SIGNAL_IDX,
  network_signal_idx = NETWORK_SIGNAL_IDX,
  network_block = NETWORK_BLOCK,
  covariate_dgp = PILOT_SETTINGS_SIGNATURE$dgp$covariates,
  network_dgp = PILOT_SETTINGS_SIGNATURE$dgp$network,
  methods = MODEL_REGISTRY,
  production_seed_table = production_seed_table,
  pilot_seed_table = pilot_seed_table,
  use_saved_competing_geometries =
    USE_SAVED_COMPETING_GEOMETRIES,
  n_loaded_complete = n_loaded_complete,
  n_computed_now = n_computed_now,
  total_production_runtime_sec =
    total_production_runtime_sec,
  metric_definitions = list(
    Overall_AUC =
      "AUC on all 400 observations against TRUE_T.",
    Hidden_AUC = paste0(
      "AUC among the 370 Y=0 observations, containing ",
      "15 hidden positives and 355 true zeros."
    ),
    Overall_LogLoss =
      "LogLoss on all 400 observations against TRUE_T.",
    Hidden_LogLoss =
      "LogLoss on the 15 hidden positives, all with TRUE_T=1.",
    Overall_Brier =
      "Brier score on all 400 observations against TRUE_T.",
    Hidden_Brier =
      "Brier score on the 15 hidden positives, all with TRUE_T=1.",
    Overall_Recall_at_45 = paste0(
      "Fraction of all 45 true positives in the top 45 ",
      "overall probability ranks."
    ),
    Hidden_Recall_at_15 = paste0(
      "Fraction of the 15 hidden positives in the top 15 ",
      "ranks among the 370 Y=0 observations."
    ),
    predictive_label_intervals = paste0(
      "Empirical equal-tail intervals from binary posterior predictive ",
      "draws generated from Bernoulli(plogis(f_i)) using quantile type 1; ",
      "compared with known TRUE_T."
    ),
    predictive_label_uncertainty = paste0(
      "For every model, unit and retained posterior beta draw, ",
      "sample T_i from Bernoulli(plogis(f_i)), with f_i = Z_i beta."
    ),
    H_interval = "not computed",
    poisson_binomial = "not used"
  ),
  predictive_interval_levels =
    PREDICTIVE_INTERVAL_LEVELS,
  predictive_quantile_type =
    PREDICTIVE_QUANTILE_TYPE
)

saveRDS(run_config, file.path(OUT_DIR, "Bayesian_linear_run_config.rds"))

save(
  global_pilot_calibration,
  production_results,
  metrics_by_replicate,
  scores_all_nodes,
  geometry_summary,
  modularity_summary,
  failure_log,
  metric_summary_long,
  main_table_means,
  paper_table_mean_sd,
  uncertainty_table_mean_sd,
  additional_diagnostic_summary,
  MODEL_REGISTRY,
  production_seed_table,
  pilot_seed_table,
  TRUE_T,
  OBSERVED_Y,
  SIGNAL_TYPE,
  TRUE_ZERO_IDX,
  HIDDEN_COV_ONLY_IDX,
  HIDDEN_NET_ONLY_IDX,
  HIDDEN_BOTH_IDX,
  OBSERVED_COV_ONLY_IDX,
  OBSERVED_NET_ONLY_IDX,
  OBSERVED_BOTH_IDX,
  HIDDEN_POS_IDX,
  OBSERVED_POS_IDX,
  TRUE_POS_IDX,
  COVARIATE_SIGNAL_IDX,
  NETWORK_SIGNAL_IDX,
  NETWORK_BLOCK,
  N_COVARIATES,
  COVARIATE_BACKGROUND_MEAN,
  COVARIATE_SIGNAL_MEAN,
  COMMON_SD,
  P_00,
  P_01,
  P_11,
  TRIAD_CLOSURE_PROB,
  run_config,
  PILOT_SETTINGS_SIGNATURE,
  PRODUCTION_SETTINGS_SIGNATURE,
  file = file.path(
    OUT_DIR,
    "SIMULATION_100_BAYESIAN_LINEAR_LINEAR_GAUSSIAN_SBM_SCAR_MIXED_SIGNAL_COMPLETE_RESULTS.RData"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo_Bayesian_linear_linear_gaussian_SBM_SCAR_mixed_signal.txt")
)

## ============================================================================
## 19. FINAL CONSOLE OUTPUT
## ============================================================================

cat("\n================ RUN ACCOUNTING ================\n")
cat("Pilot datasets          :", N_PILOT, "\n")
cat("Production replications :", N_REP, "\n")
cat("Computed in this run    :", n_computed_now, "\n")
cat("Loaded from checkpoint  :", n_loaded_complete, "\n")
cat(
  "Production runtime      :",
  format_duration(total_production_runtime_sec),
  "\n"
)

cat("\n================ FINAL MONTE CARLO MEANS ================\n")
print(
  main_table_means,
  row.names = FALSE
)

cat("\n================ FINAL PAPER TABLE ================\n")
print(
  paper_table_mean_sd,
  row.names = FALSE
)

cat("\n================ UNCERTAINTY METRICS: 80% AND 95% ================\n")
print(
  uncertainty_table_mean_sd,
  row.names = FALSE
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




## ============================================================================
## ADD SD IN PARENTHESES TO ALL COVERAGE COLUMNS
## ----------------------------------------------------------------------------
## Runs after metrics_by_replicate, method_order and paper_table_mean_sd
## have been created. No model is refitted.
## ============================================================================

required_objects_patch_A <- c(
  "metrics_by_replicate",
  "method_order",
  "paper_table_mean_sd",
  "OUT_DIR"
)

missing_objects_patch_A <- required_objects_patch_A[
  !vapply(required_objects_patch_A, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects_patch_A) > 0L) {
  stop(
    "PATCH A is missing the following objects: ",
    paste(missing_objects_patch_A, collapse = ", "),
    ". Load the completed results or run the script through Section 17 first."
  )
}

coverage_metric_to_column <- c(
  Overall_Coverage_80 = "Overall Coverage 80%",
  Hidden_Coverage_80  = "Hidden Coverage 80%",
  Overall_Coverage_95 = "Overall Coverage 95%",
  Hidden_Coverage_95  = "Hidden Coverage 95%"
)

format_coverage_mean_sd <- function(values, digits = 3L) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  
  if (length(values) == 0L) {
    return("")
  }
  
  paste0(
    formatC(mean(values), format = "f", digits = digits),
    " (",
    formatC(stats::sd(values), format = "f", digits = digits),
    ")"
  )
}

for (metric_name in names(coverage_metric_to_column)) {
  output_column <- unname(coverage_metric_to_column[metric_name])
  
  if (!metric_name %in% names(metrics_by_replicate)) {
    stop("PATCH A could not find metric: ", metric_name)
  }
  
  paper_table_mean_sd[[output_column]] <- vapply(
    method_order,
    function(method_name) {
      values_now <- metrics_by_replicate[
        metrics_by_replicate$method == method_name,
        metric_name
      ]
      
      format_coverage_mean_sd(
        values = values_now,
        digits = 3L
      )
    },
    character(1)
  )
}

write.csv(
  paper_table_mean_sd,
  file.path(
    OUT_DIR,
    "Bayesian_linear_FINAL_PAPER_TABLE_mean_SD.csv"
  ),
  row.names = FALSE,
  na = ""
)

uncertainty_columns_patch_A <- c(
  "Method",
  "Overall Coverage 80%",
  "Hidden Coverage 80%",
  "Overall Interval Score 80%",
  "Hidden Interval Score 80%",
  "Overall Coverage 95%",
  "Hidden Coverage 95%",
  "Overall Interval Score 95%",
  "Hidden Interval Score 95%"
)

if (all(uncertainty_columns_patch_A %in% names(paper_table_mean_sd))) {
  uncertainty_table_mean_sd <- paper_table_mean_sd[
    ,
    uncertainty_columns_patch_A,
    drop = FALSE
  ]
  
  write.csv(
    uncertainty_table_mean_sd,
    file.path(
      OUT_DIR,
      "Bayesian_linear_uncertainty_metrics_80_95.csv"
    ),
    row.names = FALSE,
    na = ""
  )
}

cat(
  "\n================ PAPER TABLE WITH COVERAGE MEAN (SD) ================\n"
)

print(
  paper_table_mean_sd,
  row.names = FALSE
)