## ============================================================================
## SIMULATION STUDY: COMPETING METHODS, LINEAR-GAUSSIAN + SBM, MIXED SIGNAL
## ----------------------------------------------------------------------------
## This is the linear-SCAR counterpart of the nonlinear mixed-signal simulation.
## The latent-label design, hidden/observed split and allocation of covariate,
## network and joint signals are unchanged. Only the forms of the
## covariate and network data-generating mechanisms are replaced.
##
## Fixed mixed-signal design in every replication
## ------------------------------------------------
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
## Therefore:
##   - 45 true positives in total;
##   - 30 observed positives;
##   - 15 hidden positives;
##   - 355 true zeros;
##   - 30 covariate-signal nodes;
##   - 30 network-signal nodes;
##   - 15 both-signal nodes.
##
## Linear Gaussian covariate mechanism
## -----------------------------------
## For the 30 covariate-signal nodes, namely the covariate-only and both-signal
## positives,
##
##   X_ij ~ N(0.5,1),  j=1,...,10.
##
## For all remaining 370 nodes, including the 15 network-only positives,
##
##   X_ij ~ N(0,1),  j=1,...,10.
##
## Thus network-only positives contain no covariate signal, as in the
## nonlinear mixed-signal simulation.
##
## Network mechanism
## -----------------
## The 30 network-signal nodes, namely the network-only and both-signal
## positives, form SBM block 1. All remaining 370 nodes, including the
## 15 covariate-only positives, form SBM block 0:
##
##   P(A_ij=1 | block 0 - block 0) = 0.10,
##   P(A_ij=1 | block 0 - block 1) = 0.08,
##   P(A_ij=1 | block 1 - block 1) = 0.25.
##
## There is no Chinese-restaurant layer, auxiliary local-community layer or
## triadic closure.
##
## SCAR masking
## ------------
## Each of the three signal types contains 15 positives: 5 hidden and 10
## observed. Hidden status is not used when generating X or A within a signal
## type. Hence the hidden rate is the same, 1/3, for covariate-only,
## network-only and both-signal positives.
##
## Competing methods and evaluation
## --------------------------------
## The seven competing methods, fitting procedures, point metrics, predictive
## intervals, checkpoint logic and Monte Carlo summaries are unchanged from the
## nonlinear SCAR script.
## ============================================================================

options(stringsAsFactors = FALSE)

## ============================================================================
## 0. PACKAGE CHECKS
## ============================================================================

required_namespaces <- c(
  "Matrix",
  "igraph",
  "RSpectra",
  "AdaSampling",
  "e1071",
  "INLA"
)

missing_namespaces <- required_namespaces[
  !vapply(
    required_namespaces,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
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
  library(AdaSampling)
  library(e1071)
})

## Use a single INLA thread to make numerical execution as reproducible as
## possible. Failure to set the option is non-fatal on older installations.
try(
  INLA::inla.setOption(num.threads = "1:1"),
  silent = TRUE
)

## ============================================================================
## 1. USER SETTINGS
## ============================================================================

SCRIPT_VERSION <- paste0(
  "simulation_competitors_linear_gaussian_SBM_SCAR_",
  "mixed_signal_all_metrics_v2_0"
)

MASTER_SEED <- 20260718L
N_REP <- 100L

N_NODES <- 400L

## --------------------------------------------------------------------------
## 1A. Fixed latent truth, mixed signals and SCAR observed labels
## --------------------------------------------------------------------------

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

N_COVARIATE_SIGNAL <- length(
  COVARIATE_SIGNAL_IDX
)

N_NETWORK_SIGNAL <- length(
  NETWORK_SIGNAL_IDX
)

N_BOTH_SIGNAL <- length(
  intersect(
    COVARIATE_SIGNAL_IDX,
    NETWORK_SIGNAL_IDX
  )
)

N_NETWORK_BACKGROUND <-
  N_NODES -
  N_NETWORK_SIGNAL

TRUE_T <- integer(N_NODES)
TRUE_T[TRUE_POS_IDX] <- 1L

OBSERVED_Y <- integer(N_NODES)
OBSERVED_Y[OBSERVED_POS_IDX] <- 1L

NETWORK_BLOCK <- integer(N_NODES)
NETWORK_BLOCK[NETWORK_SIGNAL_IDX] <- 1L

SIGNAL_TYPE <- rep(
  "true zero",
  N_NODES
)

SIGNAL_TYPE[HIDDEN_COV_ONLY_IDX] <-
  "hidden: covariates only"

SIGNAL_TYPE[HIDDEN_NET_ONLY_IDX] <-
  "hidden: network only"

SIGNAL_TYPE[HIDDEN_BOTH_IDX] <-
  "hidden: both"

SIGNAL_TYPE[OBSERVED_COV_ONLY_IDX] <-
  "observed: covariates only"

SIGNAL_TYPE[OBSERVED_NET_ONLY_IDX] <-
  "observed: network only"

SIGNAL_TYPE[OBSERVED_BOTH_IDX] <-
  "observed: both"

SCENARIO_NAME <- paste0(
  "Linear Gaussian covariates and Bernoulli SBM, ",
  "fixed mixed-signal SCAR design"
)

OUT_DIR <- file.path(
  getwd(),
  paste0(
    "SIMULATION_100_COMPETING_METHODS_",
    "LINEAR_GAUSSIAN_SBM_SCAR_MIXED_SIGNAL"
  )
)

CHECKPOINT_DIR <- file.path(
  OUT_DIR,
  "replicate_checkpoints"
)

GEOMETRY_DIR <- file.path(
  OUT_DIR,
  "generated_geometries"
)

dir.create(
  OUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  CHECKPOINT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  GEOMETRY_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

RESUME_REPLICATES <- TRUE
FORCE_FRESH_RUN <- FALSE
FAIL_IF_ANY_METHOD_FAILS <- TRUE

SAVE_GENERATED_GEOMETRIES <- TRUE

## --------------------------------------------------------------------------
## 1B. Evaluation settings
## --------------------------------------------------------------------------

HIDDEN_TOP_K <- 15L
OVERALL_RECALL_K <- 45L

PREDICTIVE_INTERVAL_LEVELS <- c(
  0.80,
  0.95
)

N_PREDICTIVE_DRAWS <- 2000L

PROBABILITY_EPS <- 1e-8

## --------------------------------------------------------------------------
## 1C. Linear Gaussian covariate DGP
## --------------------------------------------------------------------------

N_COVARIATES <- 10L
COVARIATE_BACKGROUND_MEAN <- 0.00
COVARIATE_SIGNAL_MEAN <- 0.50
COMMON_SD <- 1.00

## Backward-compatible aliases used in saved metadata.
ZERO_MEAN <- COVARIATE_BACKGROUND_MEAN
POSITIVE_MEAN <- COVARIATE_SIGNAL_MEAN

## --------------------------------------------------------------------------
## 1D. Two-block Bernoulli SBM
## --------------------------------------------------------------------------

P_00 <- 0.10
P_01 <- 0.08
P_11 <- 0.25

TRIAD_CLOSURE_PROB <- 0.00
TRIAD_MAX_NEW_EDGES <- Inf

## --------------------------------------------------------------------------
## 1E. Label-blind network modularity features
## --------------------------------------------------------------------------

MOD_EIG_TOL <- 1e-8
MOD_MAX_Q <- Inf
MOD_K_MAX <- 20L

## --------------------------------------------------------------------------
## 1F. Competing-method settings
## --------------------------------------------------------------------------

ADA_S <- 1
ADA_C <- 50
ADA_SAMPLE_FACTOR <- 1

WLR_ZERO_WEIGHT <- 0.5

GNN_EPOCHS <- 200L
GNN_HIDDEN <- 64L
GNN_LR <- 1e-2
GNN_WEIGHT_DECAY <- 5e-4
GNN_DROPOUT <- 0.5

## --------------------------------------------------------------------------
## 1G. Reproducibility settings
## --------------------------------------------------------------------------

RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

set.seed(MASTER_SEED)

seed_table <- data.frame(
  replicate = seq_len(N_REP),

  covariate_seed = as.integer(
    MASTER_SEED +
      100000L +
      1000L * seq_len(N_REP) +
      11L
  ),

  network_seed = as.integer(
    MASTER_SEED +
      200000L +
      1000L * seq_len(N_REP) +
      29L
  ),

  modularity_seed = as.integer(
    MASTER_SEED +
      300000L +
      1000L * seq_len(N_REP) +
      43L
  ),

  ada_seed = as.integer(
    MASTER_SEED +
      400000L +
      1000L * seq_len(N_REP) +
      101L
  ),

  wlr_seed = as.integer(
    MASTER_SEED +
      500000L +
      1000L * seq_len(N_REP) +
      103L
  ),

  inla_seed = as.integer(
    MASTER_SEED +
      600000L +
      1000L * seq_len(N_REP) +
      107L
  ),

  gnn_seed = as.integer(
    MASTER_SEED +
      700000L +
      1000L * seq_len(N_REP) +
      109L
  ),

  wlr_predictive_seed = as.integer(
    MASTER_SEED +
      800000L +
      1000L * seq_len(N_REP) +
      113L
  ),

  wlr_net_predictive_seed = as.integer(
    MASTER_SEED +
      900000L +
      1000L * seq_len(N_REP) +
      127L
  ),

  inla_predictive_seed = as.integer(
    MASTER_SEED +
      1000000L +
      1000L * seq_len(N_REP) +
      131L
  ),

  inla_net_predictive_seed = as.integer(
    MASTER_SEED +
      1100000L +
      1000L * seq_len(N_REP) +
      137L
  ),

  stringsAsFactors = FALSE
)

seed_columns <- setdiff(
  names(seed_table),
  "replicate"
)

if (
  any(
    as.matrix(
      seed_table[
        ,
        seed_columns,
        drop = FALSE
      ]
    ) >
      .Machine$integer.max
  )
) {
  stop(
    "At least one deterministic seed exceeds the R integer range."
  )
}

write.csv(
  seed_table,
  file.path(
    OUT_DIR,
    "simulation_seed_table.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 2. DESIGN CHECKS
## ============================================================================

if (N_REP != 100L) {
  stop(
    "The requested simulation study requires N_REP = 100."
  )
}

if (N_NODES != 400L) {
  stop(
    "The requested simulation study requires exactly 400 observations."
  )
}

if (!identical(TRUE_ZERO_IDX, 1:355)) {
  stop(
    "TRUE_ZERO_IDX must be exactly 1:355."
  )
}

if (!identical(HIDDEN_COV_ONLY_IDX, 356:360)) {
  stop(
    "HIDDEN_COV_ONLY_IDX must be exactly 356:360."
  )
}

if (!identical(HIDDEN_NET_ONLY_IDX, 361:365)) {
  stop(
    "HIDDEN_NET_ONLY_IDX must be exactly 361:365."
  )
}

if (!identical(HIDDEN_BOTH_IDX, 366:370)) {
  stop(
    "HIDDEN_BOTH_IDX must be exactly 366:370."
  )
}

if (!identical(OBSERVED_COV_ONLY_IDX, 371:380)) {
  stop(
    "OBSERVED_COV_ONLY_IDX must be exactly 371:380."
  )
}

if (!identical(OBSERVED_NET_ONLY_IDX, 381:390)) {
  stop(
    "OBSERVED_NET_ONLY_IDX must be exactly 381:390."
  )
}

if (!identical(OBSERVED_BOTH_IDX, 391:400)) {
  stop(
    "OBSERVED_BOTH_IDX must be exactly 391:400."
  )
}

if (!identical(HIDDEN_POS_IDX, 356:370)) {
  stop(
    "HIDDEN_POS_IDX must be exactly 356:370."
  )
}

if (!identical(OBSERVED_POS_IDX, 371:400)) {
  stop(
    "OBSERVED_POS_IDX must be exactly 371:400."
  )
}

if (!identical(TRUE_POS_IDX, 356:400)) {
  stop(
    "TRUE_POS_IDX must be exactly 356:400."
  )
}

if (
  N_TRUE_POS != 45L ||
  N_HIDDEN != 15L ||
  N_OBSERVED_POS != 30L ||
  N_TRUE_ZERO != 355L ||
  N_UNLABELED != 370L
) {
  stop(
    "The fixed 355/15/30 label design has been altered."
  )
}

if (
  length(HIDDEN_COV_ONLY_IDX) != 5L ||
  length(HIDDEN_NET_ONLY_IDX) != 5L ||
  length(HIDDEN_BOTH_IDX) != 5L ||
  length(OBSERVED_COV_ONLY_IDX) != 10L ||
  length(OBSERVED_NET_ONLY_IDX) != 10L ||
  length(OBSERVED_BOTH_IDX) != 10L
) {
  stop(
    "The required hidden 5+5+5 and observed 10+10+10 split has changed."
  )
}

if (
  N_COVARIATE_SIGNAL != 30L ||
  N_NETWORK_SIGNAL != 30L ||
  N_BOTH_SIGNAL != 15L
) {
  stop(
    "The mixed-signal design requires 30 covariate, 30 network and 15 both nodes."
  )
}

if (
  !identical(
    intersect(
      COVARIATE_SIGNAL_IDX,
      NETWORK_SIGNAL_IDX
    ),
    c(
      HIDDEN_BOTH_IDX,
      OBSERVED_BOTH_IDX
    )
  )
) {
  stop(
    "The both-signal indices are inconsistent."
  )
}

if (
  !identical(
    which(TRUE_T == 1L),
    TRUE_POS_IDX
  ) ||
  !identical(
    which(OBSERVED_Y == 1L),
    OBSERVED_POS_IDX
  ) ||
  !identical(
    which(
      TRUE_T == 1L &
        OBSERVED_Y == 0L
    ),
    HIDDEN_POS_IDX
  )
) {
  stop(
    "TRUE_T and OBSERVED_Y are inconsistent with the fixed design."
  )
}

if (
  !identical(
    which(NETWORK_BLOCK == 1L),
    NETWORK_SIGNAL_IDX
  )
) {
  stop(
    "NETWORK_BLOCK does not match NETWORK_SIGNAL_IDX."
  )
}

if (
  N_COVARIATES != 10L ||
  !is.finite(COVARIATE_BACKGROUND_MEAN) ||
  !is.finite(COVARIATE_SIGNAL_MEAN) ||
  !is.finite(COMMON_SD) ||
  COMMON_SD <= 0
) {
  stop(
    "The linear Gaussian covariate DGP is invalid."
  )
}

if (
  COVARIATE_BACKGROUND_MEAN != 0 ||
  COVARIATE_SIGNAL_MEAN != 0.5 ||
  COMMON_SD != 1
) {
  stop(
    "The requested Gaussian distributions are N(0,1) and N(0.5,1)."
  )
}

if (
  any(
    !is.finite(
      c(
        P_00,
        P_01,
        P_11
      )
    )
  ) ||
  any(
    c(
      P_00,
      P_01,
      P_11
    ) < 0
  ) ||
  any(
    c(
      P_00,
      P_01,
      P_11
    ) > 1
  )
) {
  stop(
    "P_00, P_01 and P_11 must lie in [0,1]."
  )
}




EXPECTED_DEGREE_NETWORK_BACKGROUND <-
  (
    N_NETWORK_BACKGROUND -
      1L
  ) *
  P_00 +
  N_NETWORK_SIGNAL *
  P_01

EXPECTED_DEGREE_NETWORK_SIGNAL <-
  N_NETWORK_BACKGROUND *
  P_01 +
  (
    N_NETWORK_SIGNAL -
      1L
  ) *
  P_11

EXPECTED_EDGE_COUNTS_NETWORK_BLOCKS <- c(
  "00" = choose(
    N_NETWORK_BACKGROUND,
    2L
  ) *
    P_00,

  "01" = N_NETWORK_BACKGROUND *
    N_NETWORK_SIGNAL *
    P_01,

  "11" = choose(
    N_NETWORK_SIGNAL,
    2L
  ) *
    P_11
)

if (
  HIDDEN_TOP_K != N_HIDDEN ||
  OVERALL_RECALL_K != N_TRUE_POS
) {
  stop(
    "The requested metrics require Hidden Recall@15 and Overall Recall@45."
  )
}

if (
  length(PREDICTIVE_INTERVAL_LEVELS) != 2L ||
  !isTRUE(
    all.equal(
      PREDICTIVE_INTERVAL_LEVELS,
      c(
        0.80,
        0.95
      )
    )
  ) ||
  any(
    PREDICTIVE_INTERVAL_LEVELS <= 0
  ) ||
  any(
    PREDICTIVE_INTERVAL_LEVELS >= 1
  )
) {
  stop(
    "PREDICTIVE_INTERVAL_LEVELS must equal c(0.80,0.95)."
  )
}

if (
  length(N_PREDICTIVE_DRAWS) != 1L ||
  !is.finite(N_PREDICTIVE_DRAWS) ||
  N_PREDICTIVE_DRAWS < 1000L
) {
  stop(
    "N_PREDICTIVE_DRAWS must be one integer of at least 1000."
  )
}

N_PREDICTIVE_DRAWS <- as.integer(
  N_PREDICTIVE_DRAWS
)

if (isTRUE(FORCE_FRESH_RUN)) {
  unlink(
    CHECKPOINT_DIR,
    recursive = TRUE,
    force = TRUE
  )

  dir.create(
    CHECKPOINT_DIR,
    showWarnings = FALSE,
    recursive = TRUE
  )
}

cat(
  "\n================ LINEAR MIXED-SIGNAL SCAR DESIGN ================\n"
)

cat(
  "Replications              :",
  N_REP,
  "\n"
)

cat(
  "Observations              :",
  N_NODES,
  "\n"
)

cat(
  "True zeros / positives    :",
  N_TRUE_ZERO,
  "/",
  N_TRUE_POS,
  "\n"
)

cat(
  "Observed / hidden positive:",
  N_OBSERVED_POS,
  "/",
  N_HIDDEN,
  "\n"
)

cat(
  "Hidden cov/net/both       :",
  length(HIDDEN_COV_ONLY_IDX),
  "/",
  length(HIDDEN_NET_ONLY_IDX),
  "/",
  length(HIDDEN_BOTH_IDX),
  "\n"
)

cat(
  "Observed cov/net/both     :",
  length(OBSERVED_COV_ONLY_IDX),
  "/",
  length(OBSERVED_NET_ONLY_IDX),
  "/",
  length(OBSERVED_BOTH_IDX),
  "\n"
)

cat(
  "Covariate/network/both signal:",
  N_COVARIATE_SIGNAL,
  "/",
  N_NETWORK_SIGNAL,
  "/",
  N_BOTH_SIGNAL,
  "\n"
)

cat(
  "Covariate background     : N(",
  COVARIATE_BACKGROUND_MEAN,
  ",",
  COMMON_SD,
  "^2)\n",
  sep = ""
)

cat(
  "Covariate signal         : N(",
  COVARIATE_SIGNAL_MEAN,
  ",",
  COMMON_SD,
  "^2)\n",
  sep = ""
)

cat(
  "SBM p00/p01/p11          :",
  P_00,
  "/",
  P_01,
  "/",
  P_11,
  "\n"
)

cat(
  "Network block sizes 0/1  :",
  N_NETWORK_BACKGROUND,
  "/",
  N_NETWORK_SIGNAL,
  "\n"
)

cat(
  "Expected degree block 0/1:",
  round(
    EXPECTED_DEGREE_NETWORK_BACKGROUND,
    3
  ),
  "/",
  round(
    EXPECTED_DEGREE_NETWORK_SIGNAL,
    3
  ),
  "\n"
)

cat(
  "Triadic closure          :",
  TRIAD_CLOSURE_PROB,
  "\n"
)

cat(
  "Predictive draws/method  :",
  N_PREDICTIVE_DRAWS,
  "\n"
)

## ============================================================================
## 3. GENERAL HELPERS
## ============================================================================

timestamp_now <- function() {
  format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )
}

format_duration <- function(seconds) {
  seconds <- as.numeric(seconds)
  if (!is.finite(seconds)) return(NA_character_)
  hours <- floor(seconds / 3600)
  minutes <- floor((seconds - 3600 * hours) / 60)
  seconds_remainder <- round(
    seconds - 3600 * hours - 60 * minutes
  )
  sprintf(
    "%02d:%02d:%02d",
    hours,
    minutes,
    seconds_remainder
  )
}

symmetrize <- function(M) {
  0.5 * (M + t(M))
}

standardize_columns <- function(M) {
  M <- as.matrix(M)
  center <- colMeans(M, na.rm = TRUE)
  scale_value <- apply(M, 2, stats::sd, na.rm = TRUE)
  center[!is.finite(center)] <- 0
  scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
  scaled <- sweep(M, 2, center, "-")
  scaled <- sweep(scaled, 2, scale_value, "/")
  scaled[!is.finite(scaled)] <- 0
  list(
    scaled = scaled,
    center = center,
    scale = scale_value
  )
}

sanitize_network_adjacency <- function(A) {
  A <- as.matrix(A)
  A[!is.finite(A)] <- 0
  A <- 1L * ((A + t(A)) > 0)
  diag(A) <- 0L
  A
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

validate_probability_vector <- function(
    probabilities,
    n_expected,
    method_label
) {
  probabilities <- as.numeric(probabilities)
  if (length(probabilities) != n_expected) {
    stop(
      method_label,
      " returned ", length(probabilities),
      " values; expected ", n_expected, "."
    )
  }
  if (any(!is.finite(probabilities))) {
    stop(method_label, " returned non-finite values.")
  }
  tolerance <- 1e-8
  if (
    any(probabilities < -tolerance) ||
    any(probabilities > 1 + tolerance)
  ) {
    stop(
      method_label,
      " returned values outside [0,1]. Range: [",
      min(probabilities), ", ", max(probabilities), "]."
    )
  }
  pmin(pmax(probabilities, 0), 1)
}

validate_predictive_label_draws <- function(
    predictive_T_draws,
    n_draws_expected,
    n_nodes_expected,
    method_label
) {
  if (is.null(predictive_T_draws)) return(NULL)
  predictive_T_draws <- as.matrix(predictive_T_draws)
  if (!all(dim(predictive_T_draws) == c(
    n_draws_expected,
    n_nodes_expected
  ))) {
    stop(
      method_label,
      " returned predictive-label draws with dimensions ",
      paste(dim(predictive_T_draws), collapse = " x "),
      "; expected ", n_draws_expected, " x ", n_nodes_expected, "."
    )
  }
  if (
    any(!is.finite(predictive_T_draws)) ||
    any(!predictive_T_draws %in% c(0, 1))
  ) {
    stop(method_label, " returned invalid predictive binary-label draws.")
  }
  storage.mode(predictive_T_draws) <- "integer"
  predictive_T_draws
}

overall_logloss <- function(
    probabilities,
    truth,
    eps = PROBABILITY_EPS
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
  probabilities <- as.numeric(probabilities)
  truth <- as.numeric(truth)
  mean((probabilities - truth)^2)
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

interval_score_vector <- function(L, U, truth, alpha) {
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
    level
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
    type = 1,
    names = FALSE
  )
  upper <- apply(
    predictive_T_draws,
    2L,
    stats::quantile,
    probs = 1 - alpha / 2,
    type = 1,
    names = FALSE
  )
  lower <- as.integer(lower)
  upper <- as.integer(upper)
  if (
    any(!lower %in% c(0L, 1L)) ||
    any(!upper %in% c(0L, 1L)) ||
    any(lower > upper)
  ) {
    stop("Invalid empirical interval for binary predictive draws.")
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
    interval$lower <= truth & truth <= interval$upper
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

empty_predictive_label_metrics <- function(n_nodes = N_NODES) {
  list(
    lower = rep(NA_integer_, n_nodes),
    upper = rep(NA_integer_, n_nodes),
    covered = rep(NA_real_, n_nodes),
    score = rep(NA_real_, n_nodes),
    overall_coverage = NA_real_,
    hidden_coverage = NA_real_,
    overall_interval_score = NA_real_,
    hidden_interval_score = NA_real_
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
    q025 = as.numeric(stats::quantile(x, 0.025, type = 7)),
    median = stats::median(x),
    q975 = as.numeric(stats::quantile(x, 0.975, type = 7))
  )
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
## 7. COMPETING-METHOD FITTERS
## ============================================================================

run_ada_scores <- function(
    X,
    y,
    seed
) {
  set.seed(seed)
  
  X <- as.matrix(X)
  y <- as.integer(y)
  
  n <- nrow(X)
  
  rownames(X) <- as.character(
    seq_len(n)
  )
  
  positive_rows <- as.character(
    which(
      y == 1L
    )
  )
  
  zero_rows <- as.character(
    which(
      y == 0L
    )
  )
  
  if (
    length(positive_rows) < 1L ||
    length(zero_rows) < 1L
  ) {
    stop(
      "Ada requires at least one observed positive and one zero-labelled case."
    )
  }
  
  fit <- AdaSampling::adaSample(
    Ps = positive_rows,
    Ns = zero_rows,
    train.mat = X,
    test.mat = X,
    classifier = "svm",
    s = ADA_S,
    C = ADA_C,
    sampleFactor = ADA_SAMPLE_FACTOR
  )
  
  validate_probability_vector(
    probabilities = fit[, 1L],
    n_expected = n,
    method_label = "AdaSampling"
  )
}

run_wlr_scores <- function(
    X,
    y,
    seed,
    predictive_seed,
    zero_weight = WLR_ZERO_WEIGHT,
    n_predictive_draws = N_PREDICTIVE_DRAWS
) {
  set.seed(seed)
  
  X <- as.matrix(X)
  y <- as.integer(y)
  n <- nrow(X)
  
  data_fit <- as.data.frame(X)
  names(data_fit) <- paste0("x", seq_len(ncol(X)))
  data_fit$Y <- y
  
  fit_weights <- ifelse(y == 1L, 1, zero_weight)
  
  fit <- suppressWarnings(
    tryCatch(
      glm(
        Y ~ .,
        data = data_fit,
        family = binomial(link = "logit"),
        weights = fit_weights,
        x = TRUE,
        model = TRUE
      ),
      error = function(e) NULL
    )
  )
  
  if (is.null(fit)) stop("WLR fit failed.")
  
  probabilities <- suppressWarnings(
    as.numeric(
      predict(
        fit,
        newdata = data_fit,
        type = "response"
      )
    )
  )
  
  probabilities <- validate_probability_vector(
    probabilities = probabilities,
    n_expected = n,
    method_label = "WLR"
  )
  
  beta_hat_full <- stats::coef(fit)
  estimable <- is.finite(beta_hat_full)
  if (!any(estimable)) {
    stop("WLR returned no finite estimable coefficients.")
  }
  
  X_design_full <- stats::model.matrix(fit)
  V_full <- suppressWarnings(stats::vcov(fit, complete = TRUE))
  
  beta_hat <- as.numeric(beta_hat_full[estimable])
  X_design <- X_design_full[, estimable, drop = FALSE]
  V_beta <- as.matrix(V_full[estimable, estimable, drop = FALSE])
  
  if (
    any(!is.finite(beta_hat)) ||
    any(!is.finite(X_design)) ||
    any(!is.finite(V_beta))
  ) {
    stop("WLR asymptotic coefficient uncertainty contains non-finite values.")
  }
  
  V_beta <- symmetrize(V_beta)
  eigen_V <- eigen(V_beta, symmetric = TRUE)
  eigen_floor <- max(1e-12, 1e-10 * max(1, max(abs(eigen_V$values))))
  V_beta <- symmetrize(
    eigen_V$vectors %*%
      (pmax(eigen_V$values, eigen_floor) * t(eigen_V$vectors))
  )
  
  eta_mean <- as.vector(X_design %*% beta_hat)
  eta_variance <- rowSums((X_design %*% V_beta) * X_design)
  eta_variance <- pmax(eta_variance, 0)
  eta_sd <- sqrt(eta_variance)
  
  set.seed(as.integer(predictive_seed))
  eta_draws <- matrix(
    stats::rnorm(
      n_predictive_draws * n,
      mean = rep(eta_mean, each = n_predictive_draws),
      sd = rep(eta_sd, each = n_predictive_draws)
    ),
    nrow = n_predictive_draws,
    ncol = n
  )
  
  probability_draws <- plogis(eta_draws)
  predictive_T_draws <- matrix(
    stats::rbinom(
      n_predictive_draws * n,
      size = 1L,
      prob = as.vector(probability_draws)
    ),
    nrow = n_predictive_draws,
    ncol = n
  )
  
  list(
    probabilities = probabilities,
    predictive_T_draws = predictive_T_draws,
    uncertainty_type = paste0(
      "WLR asymptotic Gaussian linear-predictor uncertainty plus ",
      "Bernoulli posterior-predictive variation"
    )
  )
}

run_inla_scores <- function(
    X,
    y,
    seed,
    predictive_seed,
    n_predictive_draws = N_PREDICTIVE_DRAWS
) {
  set.seed(seed)
  
  X <- as.matrix(X)
  y <- as.integer(y)
  n <- nrow(X)
  d <- ncol(X)
  
  covariate_names <- paste0("x", seq_len(d))
  data_features <- as.data.frame(X)
  names(data_features) <- covariate_names
  data_features$Intercept <- 1
  
  inla_formula <- as.formula(
    paste0(
      "y ~ 0 + Intercept + ",
      paste(covariate_names, collapse = " + ")
    )
  )
  
  inla_data <- cbind(
    data.frame(
      y = y,
      Ntrials = rep(1L, n)
    ),
    data_features
  )
  
  fit <- tryCatch(
    INLA::inla(
      inla_formula,
      family = "binomial",
      data = inla_data,
      Ntrials = inla_data$Ntrials,
      control.predictor = list(
        compute = TRUE
      ),
      control.compute = list(
        dic = FALSE,
        waic = FALSE,
        cpo = FALSE,
        return.marginals.predictor = TRUE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) stop("INLA fit failed.")
  if (is.null(fit$summary.linear.predictor)) {
    stop("INLA did not return summary.linear.predictor.")
  }
  if (
    is.null(fit$marginals.linear.predictor) ||
    length(fit$marginals.linear.predictor) != n
  ) {
    stop(
      "INLA did not return one posterior linear-predictor marginal per node."
    )
  }
  
  linear_predictor_mean <- as.numeric(
    fit$summary.linear.predictor[, "mean"]
  )
  probabilities <- plogis(linear_predictor_mean)
  probabilities <- validate_probability_vector(
    probabilities = probabilities,
    n_expected = n,
    method_label = "INLA"
  )
  
  set.seed(as.integer(predictive_seed))
  predictive_T_draws <- matrix(
    NA_integer_,
    nrow = n_predictive_draws,
    ncol = n
  )
  
  for (ii in seq_len(n)) {
    eta_draws_i <- INLA::inla.rmarginal(
      n = n_predictive_draws,
      marginal = fit$marginals.linear.predictor[[ii]]
    )
    probability_draws_i <- plogis(eta_draws_i)
    predictive_T_draws[, ii] <- stats::rbinom(
      n_predictive_draws,
      size = 1L,
      prob = probability_draws_i
    )
  }
  
  list(
    probabilities = probabilities,
    predictive_T_draws = predictive_T_draws,
    uncertainty_type = paste0(
      "INLA posterior linear-predictor marginals plus Bernoulli ",
      "posterior-predictive variation"
    )
  )
}

renorm_adj <- function(A) {
  A <- sanitize_network_adjacency(
    A
  )
  
  n <- nrow(A)
  
  A_tilde <- A +
    diag(
      1,
      n
    )
  
  degree <- rowSums(
    A_tilde
  )
  
  D_inverse_sqrt <- diag(
    ifelse(
      degree > 0,
      1 / sqrt(degree),
      0
    ),
    n
  )
  
  D_inverse_sqrt %*%
    A_tilde %*%
    D_inverse_sqrt
}

run_gnn_scores <- function(
    X,
    y,
    A,
    seed,
    epochs = GNN_EPOCHS,
    hidden = GNN_HIDDEN,
    lr = GNN_LR,
    weight_decay = GNN_WEIGHT_DECAY,
    dropout = GNN_DROPOUT
) {
  set.seed(seed)
  
  X <- as.matrix(X)
  y <- as.numeric(y)
  
  n <- nrow(X)
  d <- ncol(X)
  
  A_hat <- renorm_adj(
    A
  )
  
  W1 <- matrix(
    rnorm(
      d * hidden,
      sd = 0.1
    ),
    d,
    hidden
  )
  
  W2 <- matrix(
    rnorm(
      hidden,
      sd = 0.1
    ),
    hidden,
    1L
  )
  
  beta1 <- 0.9
  beta2 <- 0.999
  eps <- 1e-8
  
  mW1 <- matrix(
    0,
    d,
    hidden
  )
  
  vW1 <- matrix(
    0,
    d,
    hidden
  )
  
  mW2 <- matrix(
    0,
    hidden,
    1L
  )
  
  vW2 <- matrix(
    0,
    hidden,
    1L
  )
  
  keep_probability <- 1 -
    dropout
  
  for (
    epoch in seq_len(epochs)
  ) {
    Y1 <- X %*%
      W1
    
    Z1 <- A_hat %*%
      Y1
    
    H <- pmax(
      Z1,
      0
    )
    
    if (dropout > 0) {
      dropout_mask <- matrix(
        rbinom(
          n * hidden,
          1,
          keep_probability
        ),
        n,
        hidden
      )
      
      H_drop <- H *
        dropout_mask /
        keep_probability
    } else {
      dropout_mask <- matrix(
        1,
        n,
        hidden
      )
      
      H_drop <- H
    }
    
    Y2 <- H_drop %*%
      W2
    
    Z2 <- A_hat %*%
      Y2
    
    probability <- plogis(
      Z2
    )[, 1L]
    
    gradient_Z2 <- matrix(
      (
        probability - y
      ) / n,
      n,
      1L
    )
    
    gradient_Y2 <- A_hat %*%
      gradient_Z2
    
    gradient_W2 <- t(
      H_drop
    ) %*%
      gradient_Y2 +
      weight_decay *
      W2
    
    gradient_H_drop <- gradient_Y2 %*%
      t(W2)
    
    gradient_H <- if (
      dropout > 0
    ) {
      gradient_H_drop *
        dropout_mask /
        keep_probability
    } else {
      gradient_H_drop
    }
    
    gradient_Z1 <- gradient_H *
      (
        Z1 > 0
      )
    
    gradient_Y1 <- A_hat %*%
      gradient_Z1
    
    gradient_W1 <- t(X) %*%
      gradient_Y1 +
      weight_decay *
      W1
    
    mW1 <- beta1 *
      mW1 +
      (
        1 - beta1
      ) *
      gradient_W1
    
    vW1 <- beta2 *
      vW1 +
      (
        1 - beta2
      ) *
      (
        gradient_W1 *
          gradient_W1
      )
    
    mW2 <- beta1 *
      mW2 +
      (
        1 - beta1
      ) *
      gradient_W2
    
    vW2 <- beta2 *
      vW2 +
      (
        1 - beta2
      ) *
      (
        gradient_W2 *
          gradient_W2
      )
    
    mW1_hat <- mW1 /
      (
        1 -
          beta1^epoch
      )
    
    vW1_hat <- vW1 /
      (
        1 -
          beta2^epoch
      )
    
    mW2_hat <- mW2 /
      (
        1 -
          beta1^epoch
      )
    
    vW2_hat <- vW2 /
      (
        1 -
          beta2^epoch
      )
    
    W1 <- W1 -
      lr *
      mW1_hat /
      (
        sqrt(vW1_hat) +
          eps
      )
    
    W2 <- W2 -
      lr *
      mW2_hat /
      (
        sqrt(vW2_hat) +
          eps
      )
  }
  
  Y1 <- X %*%
    W1
  
  Z1 <- A_hat %*%
    Y1
  
  H <- pmax(
    Z1,
    0
  )
  
  Y2 <- H %*%
    W2
  
  Z2 <- A_hat %*%
    Y2
  
  probabilities <- plogis(
    Z2
  )[, 1L]
  
  validate_probability_vector(
    probabilities = probabilities,
    n_expected = n,
    method_label = "GNN"
  )
}

method_registry <- data.frame(
  method_id = c(
    "ADA",
    "ADA_NET",
    "WLR",
    "WLR_NET",
    "INLA",
    "INLA_NET",
    "GNN"
  ),
  method = c(
    "Ada",
    "Ada + Net",
    "WLR",
    "WLR + Net",
    "INLA",
    "INLA + Net",
    "GNN"
  ),
  stringsAsFactors = FALSE
)

run_competing_method <- function(
    method_id,
    X_cov,
    X_cov_net,
    A_network,
    y,
    replicate_seed_row
) {
  switch(
    method_id,
    
    ADA = list(
      probabilities = run_ada_scores(
        X = X_cov,
        y = y,
        seed = replicate_seed_row$ada_seed
      ),
      predictive_T_draws = NULL,
      uncertainty_type = "not available"
    ),
    
    ADA_NET = list(
      probabilities = run_ada_scores(
        X = X_cov_net,
        y = y,
        seed = replicate_seed_row$ada_seed
      ),
      predictive_T_draws = NULL,
      uncertainty_type = "not available"
    ),
    
    WLR = run_wlr_scores(
      X = X_cov,
      y = y,
      seed = replicate_seed_row$wlr_seed,
      predictive_seed = replicate_seed_row$wlr_predictive_seed
    ),
    
    WLR_NET = run_wlr_scores(
      X = X_cov_net,
      y = y,
      seed = replicate_seed_row$wlr_seed,
      predictive_seed = replicate_seed_row$wlr_net_predictive_seed
    ),
    
    INLA = run_inla_scores(
      X = X_cov,
      y = y,
      seed = replicate_seed_row$inla_seed,
      predictive_seed = replicate_seed_row$inla_predictive_seed
    ),
    
    INLA_NET = run_inla_scores(
      X = X_cov_net,
      y = y,
      seed = replicate_seed_row$inla_seed,
      predictive_seed = replicate_seed_row$inla_net_predictive_seed
    ),
    
    GNN = list(
      probabilities = run_gnn_scores(
        X = X_cov,
        y = y,
        A = A_network,
        seed = replicate_seed_row$gnn_seed
      ),
      predictive_T_draws = NULL,
      uncertainty_type = "not available"
    ),
    
    stop("Unknown method_id: ", method_id)
  )
}


## ============================================================================
## 8. COMMON EVALUATION
## ============================================================================

evaluate_method_scores <- function(
    method_output,
    method_id,
    method_label,
    replicate_number
) {
  if (!is.list(method_output) || is.null(method_output$probabilities)) {
    stop(method_label, " did not return the required method-output list.")
  }
  
  probabilities <- validate_probability_vector(
    probabilities = method_output$probabilities,
    n_expected = N_NODES,
    method_label = method_label
  )
  
  predictive_T_draws <- validate_predictive_label_draws(
    predictive_T_draws = method_output$predictive_T_draws,
    n_draws_expected = N_PREDICTIVE_DRAWS,
    n_nodes_expected = N_NODES,
    method_label = method_label
  )
  
  uncertainty_type <- if (
    is.null(method_output$uncertainty_type) ||
    length(method_output$uncertainty_type) != 1L
  ) {
    NA_character_
  } else {
    as.character(method_output$uncertainty_type)
  }
  
  unlabeled_idx <- which(OBSERVED_Y == 0L)
  hidden_truth_unlabeled <- as.integer(
    TRUE_T[unlabeled_idx] == 1L
  )
  
  if (
    length(unlabeled_idx) != N_UNLABELED ||
    sum(hidden_truth_unlabeled) != N_HIDDEN
  ) {
    stop("The unlabeled evaluation set is inconsistent.")
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
  
  if (is.null(predictive_T_draws)) {
    predictive_80 <- empty_predictive_label_metrics(N_NODES)
    predictive_95 <- empty_predictive_label_metrics(N_NODES)
  } else {
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
  }
  
  metrics <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    method_id = method_id,
    method = method_label,
    uncertainty_type = uncertainty_type,
    
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
    Overall_Interval_Score_80 = predictive_80$overall_interval_score,
    Hidden_Interval_Score_80 = predictive_80$hidden_interval_score,
    
    Overall_Coverage_95 = predictive_95$overall_coverage,
    Hidden_Coverage_95 = predictive_95$hidden_coverage,
    Overall_Interval_Score_95 = predictive_95$overall_interval_score,
    Hidden_Interval_Score_95 = predictive_95$hidden_interval_score,
    
    stringsAsFactors = FALSE
  )
  
  node_scores <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    method_id = method_id,
    method = method_label,
    uncertainty_type = uncertainty_type,
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
## 9. SETTINGS SIGNATURE FOR CHECKPOINT SAFETY
## ============================================================================

SETTINGS_SIGNATURE <- list(
  script_version = SCRIPT_VERSION,
  master_seed = MASTER_SEED,
  n_rep = N_REP,
  n_nodes = N_NODES,

  latent_truth_and_signal_allocation = list(
    true_zero_idx = TRUE_ZERO_IDX,

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

    hidden_positive_idx =
      HIDDEN_POS_IDX,

    observed_positive_idx =
      OBSERVED_POS_IDX,

    true_positive_idx =
      TRUE_POS_IDX,

    covariate_signal_idx =
      COVARIATE_SIGNAL_IDX,

    network_signal_idx =
      NETWORK_SIGNAL_IDX
  ),

  label_mechanism = list(
    type = "SCAR with fixed balanced signal-type allocation",

    construction = paste0(
      "Each of covariate-only, network-only and both-signal groups ",
      "contains 15 positives, with 5 hidden and 10 observed. ",
      "Hidden status is not used in X or A generation within signal type."
    ),

    hidden_fraction_within_each_signal_type =
      1 / 3,

    n_hidden =
      N_HIDDEN,

    n_observed_positive =
      N_OBSERVED_POS
  ),

  covariate_dgp = list(
    type = paste0(
      "ten independent Gaussian covariates determined by ",
      "covariate-signal membership"
    ),

    n_covariates =
      N_COVARIATES,

    background_mean =
      COVARIATE_BACKGROUND_MEAN,

    signal_mean =
      COVARIATE_SIGNAL_MEAN,

    common_sd =
      COMMON_SD,

    background_contains = c(
      "355 true zeros",
      "15 network-only positives"
    ),

    signal_contains = c(
      "15 covariate-only positives",
      "15 both-signal positives"
    )
  ),

  network_dgp = list(
    type = paste0(
      "two-block Bernoulli SBM determined by ",
      "network-signal membership"
    ),

    block_zero_size =
      N_NETWORK_BACKGROUND,

    block_one_size =
      N_NETWORK_SIGNAL,

    block_zero_contains = c(
      "355 true zeros",
      "15 covariate-only positives"
    ),

    block_one_contains = c(
      "15 network-only positives",
      "15 both-signal positives"
    ),

    p_00 = P_00,
    p_01 = P_01,
    p_11 = P_11,

    triad_closure_prob =
      TRIAD_CLOSURE_PROB,

    triad_max_new_edges =
      TRIAD_MAX_NEW_EDGES,

    expected_degree_network_background =
      EXPECTED_DEGREE_NETWORK_BACKGROUND,

    expected_degree_network_signal =
      EXPECTED_DEGREE_NETWORK_SIGNAL
  ),

  modularity = list(
    eig_tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q,
    k_max = MOD_K_MAX
  ),

  methods = list(
    registry = method_registry,
    ada_s = ADA_S,
    ada_C = ADA_C,
    ada_sample_factor = ADA_SAMPLE_FACTOR,
    wlr_zero_weight = WLR_ZERO_WEIGHT,
    gnn_epochs = GNN_EPOCHS,
    gnn_hidden = GNN_HIDDEN,
    gnn_lr = GNN_LR,
    gnn_weight_decay = GNN_WEIGHT_DECAY,
    gnn_dropout = GNN_DROPOUT
  ),

  evaluation = list(
    hidden_top_k = HIDDEN_TOP_K,
    overall_recall_k = OVERALL_RECALL_K,
    predictive_interval_levels =
      PREDICTIVE_INTERVAL_LEVELS,
    n_predictive_draws =
      N_PREDICTIVE_DRAWS,
    interval_quantile_type = 1L,

    truth_for_evaluation = paste0(
      "known latent TRUE_T; hidden positives count as T=1"
    ),

    uncertainty_available_for = c(
      "WLR",
      "WLR + Net",
      "INLA",
      "INLA + Net"
    ),

    uncertainty_not_available_for = c(
      "Ada",
      "Ada + Net",
      "GNN"
    ),

    H_interval = "not computed",
    poisson_binomial = "not used"
  )
)

## ============================================================================
## 10. ONE COMPLETE SIMULATION REPLICATION
## ============================================================================

run_one_replication <- function(
    replicate_number
) {
  replicate_number <- as.integer(
    replicate_number
  )
  
  seed_row <- seed_table[
    seed_table$replicate ==
      replicate_number,
    ,
    drop = FALSE
  ]
  
  if (nrow(seed_row) != 1L) {
    stop(
      "Seed-table lookup failed for replicate ",
      replicate_number,
      "."
    )
  }
  
  checkpoint_file <- file.path(
    CHECKPOINT_DIR,
    sprintf(
      "replicate_%03d.rds",
      replicate_number
    )
  )
  
  checkpoint <- NULL
  
  if (
    isTRUE(RESUME_REPLICATES) &&
    !isTRUE(FORCE_FRESH_RUN) &&
    file.exists(checkpoint_file)
  ) {
    checkpoint_candidate <- tryCatch(
      readRDS(
        checkpoint_file
      ),
      error = function(e) {
        NULL
      }
    )
    
    if (
      is.list(checkpoint_candidate) &&
      identical(
        checkpoint_candidate$settings_signature,
        SETTINGS_SIGNATURE
      ) &&
      identical(
        as.integer(
          checkpoint_candidate$replicate
        ),
        replicate_number
      )
    ) {
      checkpoint <- checkpoint_candidate
    }
  }
  
  if (
    !is.null(checkpoint) &&
    isTRUE(
      checkpoint$complete
    )
  ) {
    return(
      checkpoint
    )
  }
  
  replication_start <- Sys.time()
  
  ## All mixed-signal linear-Gaussian/SBM geometry is regenerated from deterministic seeds when a partial
  ## checkpoint is resumed. This guarantees that every method sees the same
  ## data without storing large temporary objects inside the checkpoint.
  covariate_set <- make_covariate_set(
    seed = seed_row$covariate_seed,
    signal_idx = COVARIATE_SIGNAL_IDX
  )
  
  network_set <- generate_network_set(
    seed = seed_row$network_seed,
    signal_idx = NETWORK_SIGNAL_IDX
  )
  
  if (
    !identical(
      sort(
        covariate_set$signal_idx
      ),
      sort(
        COVARIATE_SIGNAL_IDX
      )
    )
  ) {
    stop(
      "Covariate-signal positions changed unexpectedly."
    )
  }
  
  if (
    !identical(
      sort(
        network_set$signal_idx
      ),
      sort(
        NETWORK_SIGNAL_IDX
      )
    )
  ) {
    stop(
      "Network-signal positions changed unexpectedly."
    )
  }
  
  modularity <- modularity_features_single_network(
    A = network_set$A,
    seed = seed_row$modularity_seed,
    eig_tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q,
    k_max = MOD_K_MAX
  )
  
  X_cov <- as.matrix(
    covariate_set$X_std
  )
  
  Z_net <- as.matrix(
    modularity$Z
  )
  
  X_cov_net <- cbind(
    X_cov,
    Z_net
  )
  
  if (
    nrow(X_cov) != N_NODES ||
    nrow(Z_net) != N_NODES ||
    nrow(X_cov_net) != N_NODES
  ) {
    stop(
      "Feature matrices do not have 400 rows."
    )
  }
  
  if (
    any(!is.finite(X_cov)) ||
    any(!is.finite(Z_net)) ||
    any(!is.finite(X_cov_net))
  ) {
    stop(
      "Non-finite values were found in the generated features."
    )
  }
  
  if (is.null(checkpoint)) {
    checkpoint <- list(
      settings_signature = SETTINGS_SIGNATURE,
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
  
  geometry_summary <- cbind(
    data.frame(
      scenario = SCENARIO_NAME,
      replicate = replicate_number,
      covariate_seed = seed_row$covariate_seed,
      network_seed = seed_row$network_seed,
      modularity_seed = seed_row$modularity_seed,
      n = N_NODES,
      n_true_zero = N_TRUE_ZERO,
      n_true_positive = N_TRUE_POS,
      n_observed_positive = N_OBSERVED_POS,
      n_hidden_positive = N_HIDDEN,
      n_covariate_signal = length(
        COVARIATE_SIGNAL_IDX
      ),
      n_network_signal = length(
        NETWORK_SIGNAL_IDX
      ),
      n_both_signal = length(
        intersect(
          COVARIATE_SIGNAL_IDX,
          NETWORK_SIGNAL_IDX
        )
      ),
      stringsAsFactors = FALSE
    ),
    network_set$summary
  )
  
  modularity_summary <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    q_network = modularity$q,
    selected_eigenvalues = paste(
      signif(
        modularity$selected_eigenvalues,
        8
      ),
      collapse = ";"
    ),
    n_covariates = ncol(X_cov),
    n_network_features = ncol(Z_net),
    n_covariates_plus_network = ncol(X_cov_net),
    stringsAsFactors = FALSE
  )
  
  checkpoint$geometry_summary <- geometry_summary
  checkpoint$modularity_summary <- modularity_summary
  
  if (isTRUE(SAVE_GENERATED_GEOMETRIES)) {
    geometry_file <- file.path(
      GEOMETRY_DIR,
      sprintf(
        "replicate_%03d_geometry.rds",
        replicate_number
      )
    )
    
    saveRDS(
      list(
        settings_signature = SETTINGS_SIGNATURE,
        scenario = SCENARIO_NAME,
        replicate = replicate_number,
        seeds = seed_row,
        T = TRUE_T,
        Y = OBSERVED_Y,
        hidden_positive_idx = HIDDEN_POS_IDX,
        observed_positive_idx = OBSERVED_POS_IDX,
        hidden_covariate_only_idx = HIDDEN_COV_ONLY_IDX,
        hidden_network_only_idx = HIDDEN_NET_ONLY_IDX,
        hidden_both_idx = HIDDEN_BOTH_IDX,
        observed_covariate_only_idx = OBSERVED_COV_ONLY_IDX,
        observed_network_only_idx = OBSERVED_NET_ONLY_IDX,
        observed_both_idx = OBSERVED_BOTH_IDX,
        signal_type = SIGNAL_TYPE,
        covariate_signal_idx = COVARIATE_SIGNAL_IDX,
        network_signal_idx = NETWORK_SIGNAL_IDX,
        network_block = NETWORK_BLOCK,
        covariate_dgp = list(
          n_covariates = N_COVARIATES,
          background_mean = COVARIATE_BACKGROUND_MEAN,
          signal_mean = COVARIATE_SIGNAL_MEAN,
          common_sd = COMMON_SD,
          allocation = "covariate-only plus both are signal"
        ),
        network_dgp = list(
          type = "two-block Bernoulli SBM",
          block_definition = "network-only plus both form block 1",
          p_00 = P_00,
          p_01 = P_01,
          p_11 = P_11,
          triad_closure_prob = TRIAD_CLOSURE_PROB
        ),
        covariate_profile = covariate_set$covariate_profile,
        X_raw = covariate_set$X_raw,
        X_std = X_cov,
        X_center = covariate_set$center,
        X_scale = covariate_set$scale,
        A_network = network_set$A,
        network_summary = network_set$summary,
        Z_network_modularity = Z_net,
        modularity_q = modularity$q,
        modularity_selected_eigenvalues =
          modularity$selected_eigenvalues
      ),
      geometry_file
    )
  }
  
  for (
    method_index in seq_len(
      nrow(method_registry)
    )
  ) {
    method_id <- method_registry$method_id[
      method_index
    ]
    
    method_label <- method_registry$method[
      method_index
    ]
    
    if (
      method_id %in%
      names(
        checkpoint$method_results
      )
    ) {
      cat(
        sprintf(
          "[%s] Replicate %03d/%03d | reusing %s\n",
          timestamp_now(),
          replicate_number,
          N_REP,
          method_label
        )
      )
      
      next
    }
    
    cat(
      sprintf(
        "[%s] Replicate %03d/%03d | fitting %s\n",
        timestamp_now(),
        replicate_number,
        N_REP,
        method_label
      )
    )
    
    method_start <- Sys.time()
    
    method_error <- NA_character_
    
    method_output <- tryCatch(
      run_competing_method(
        method_id = method_id,
        X_cov = X_cov,
        X_cov_net = X_cov_net,
        A_network = network_set$A,
        y = OBSERVED_Y,
        replicate_seed_row = seed_row
      ),
      error = function(e) {
        method_error <<- conditionMessage(e)
        NULL
      }
    )
    
    method_runtime_sec <- as.numeric(
      difftime(
        Sys.time(),
        method_start,
        units = "secs"
      )
    )
    
    if (is.null(method_output)) {
      checkpoint$method_failures[[method_id]] <- data.frame(
        scenario = SCENARIO_NAME,
        replicate = replicate_number,
        method_id = method_id,
        method = method_label,
        runtime_sec = method_runtime_sec,
        error_message = method_error,
        stringsAsFactors = FALSE
      )
      
      saveRDS(checkpoint, checkpoint_file)
      
      if (isTRUE(FAIL_IF_ANY_METHOD_FAILS)) {
        stop(
          "Method failure in replicate ",
          replicate_number,
          " for ",
          method_label,
          ": ",
          method_error,
          ". Partial checkpoint saved to ",
          checkpoint_file,
          "."
        )
      }
      
      next
    }
    
    evaluation <- evaluate_method_scores(
      method_output = method_output,
      method_id = method_id,
      method_label = method_label,
      replicate_number = replicate_number
    )
    
    evaluation$metrics$runtime_sec <- method_runtime_sec
    
    checkpoint$method_results[[method_id]] <- evaluation
    checkpoint$method_failures[[method_id]] <- NULL
    
    saveRDS(
      checkpoint,
      checkpoint_file
    )
    
    rm(
      method_output,
      evaluation
    )
    
    invisible(
      gc()
    )
  }
  
  completed_method_ids <- names(
    checkpoint$method_results
  )
  
  missing_method_ids <- setdiff(
    method_registry$method_id,
    completed_method_ids
  )
  
  if (length(missing_method_ids) > 0L) {
    checkpoint$complete <- FALSE
    
    saveRDS(
      checkpoint,
      checkpoint_file
    )
    
    stop(
      "Replication ",
      replicate_number,
      " did not complete every method. Missing: ",
      paste(
        missing_method_ids,
        collapse = ", "
      ),
      "."
    )
  }
  
  checkpoint$complete <- TRUE
  checkpoint$time_end <- Sys.time()
  
  checkpoint$runtime_sec <- as.numeric(
    difftime(
      checkpoint$time_end,
      checkpoint$time_start,
      units = "secs"
    )
  )
  
  saveRDS(
    checkpoint,
    checkpoint_file
  )
  
  rm(
    covariate_set,
    network_set,
    modularity,
    X_cov,
    Z_net,
    X_cov_net
  )
  
  invisible(
    gc()
  )
  
  checkpoint
}

## ============================================================================
## 11. RUN ALL 100 REPLICATIONS
## ============================================================================

replication_results <- vector(
  "list",
  N_REP
)

simulation_start <- Sys.time()

n_loaded_complete <- 0L
n_computed_now <- 0L

for (
  replicate_number in seq_len(N_REP)
) {
  checkpoint_file <- file.path(
    CHECKPOINT_DIR,
    sprintf(
      "replicate_%03d.rds",
      replicate_number
    )
  )
  
  was_complete_before <- FALSE
  
  if (
    isTRUE(RESUME_REPLICATES) &&
    !isTRUE(FORCE_FRESH_RUN) &&
    file.exists(checkpoint_file)
  ) {
    prior_checkpoint <- tryCatch(
      readRDS(
        checkpoint_file
      ),
      error = function(e) {
        NULL
      }
    )
    
    was_complete_before <-
      is.list(prior_checkpoint) &&
      identical(
        prior_checkpoint$settings_signature,
        SETTINGS_SIGNATURE
      ) &&
      isTRUE(
        prior_checkpoint$complete
      )
  }
  
  replication_results[[replicate_number]] <- run_one_replication(
    replicate_number
  )
  
  if (was_complete_before) {
    n_loaded_complete <- n_loaded_complete +
      1L
  } else {
    n_computed_now <- n_computed_now +
      1L
  }
  
  elapsed_sec <- as.numeric(
    difftime(
      Sys.time(),
      simulation_start,
      units = "secs"
    )
  )
  
  mean_sec_per_processed <- elapsed_sec /
    replicate_number
  
  projected_total_sec <- mean_sec_per_processed *
    N_REP
  
  cat(
    sprintf(
      "[%s] Completed replication %03d/%03d | elapsed %s | projected total %s\n",
      timestamp_now(),
      replicate_number,
      N_REP,
      format_duration(
        elapsed_sec
      ),
      format_duration(
        projected_total_sec
      )
    )
  )
}

simulation_end <- Sys.time()

## ============================================================================
## 12. COMBINE RESULTS
## ============================================================================

if (
  length(replication_results) != N_REP ||
  any(
    !vapply(
      replication_results,
      function(x) {
        is.list(x) &&
          isTRUE(x$complete)
      },
      logical(1)
    )
  )
) {
  stop(
    "Not all 100 simulation replications completed successfully."
  )
}

metrics_by_replicate <- do.call(
  rbind,
  lapply(
    replication_results,
    function(replication_result) {
      do.call(
        rbind,
        lapply(
          method_registry$method_id,
          function(method_id) {
            replication_result$method_results[[method_id]]$metrics
          }
        )
      )
    }
  )
)

row.names(
  metrics_by_replicate
) <- NULL

scores_all_nodes <- do.call(
  rbind,
  lapply(
    replication_results,
    function(replication_result) {
      do.call(
        rbind,
        lapply(
          method_registry$method_id,
          function(method_id) {
            replication_result$method_results[[method_id]]$node_scores
          }
        )
      )
    }
  )
)

row.names(
  scores_all_nodes
) <- NULL

geometry_summary <- do.call(
  rbind,
  lapply(
    replication_results,
    `[[`,
    "geometry_summary"
  )
)

row.names(
  geometry_summary
) <- NULL

modularity_summary <- do.call(
  rbind,
  lapply(
    replication_results,
    `[[`,
    "modularity_summary"
  )
)

row.names(
  modularity_summary
) <- NULL

failure_log_parts <- lapply(
  replication_results,
  function(replication_result) {
    failures_now <- replication_result$method_failures
    
    if (length(failures_now) == 0L) {
      return(NULL)
    }
    
    do.call(
      rbind,
      failures_now
    )
  }
)

failure_log_parts <- Filter(
  Negate(is.null),
  failure_log_parts
)

failure_log <- if (
  length(failure_log_parts) == 0L
) {
  data.frame(
    scenario = character(0),
    replicate = integer(0),
    method_id = character(0),
    method = character(0),
    runtime_sec = numeric(0),
    error_message = character(0),
    stringsAsFactors = FALSE
  )
} else {
  do.call(
    rbind,
    failure_log_parts
  )
}

if (
  nrow(metrics_by_replicate) !=
  N_REP *
  nrow(method_registry)
) {
  stop(
    "The combined metric table does not contain 100 x 7 rows."
  )
}

if (
  nrow(scores_all_nodes) !=
  N_REP *
  nrow(method_registry) *
  N_NODES
) {
  stop(
    "The combined node-score table does not contain 100 x 7 x 400 rows."
  )
}

method_order <- method_registry$method

metrics_by_replicate$method <- factor(
  metrics_by_replicate$method,
  levels = method_order
)

metrics_by_replicate <- metrics_by_replicate[
  order(
    metrics_by_replicate$method,
    metrics_by_replicate$replicate
  ),
  ,
  drop = FALSE
]

metrics_by_replicate$method <- as.character(
  metrics_by_replicate$method
)

scores_all_nodes$method <- factor(
  scores_all_nodes$method,
  levels = method_order
)

scores_all_nodes <- scores_all_nodes[
  order(
    scores_all_nodes$method,
    scores_all_nodes$replicate,
    scores_all_nodes$node
  ),
  ,
  drop = FALSE
]

scores_all_nodes$method <- as.character(
  scores_all_nodes$method
)

write.csv(
  metrics_by_replicate,
  file.path(
    OUT_DIR,
    "competing_methods_metrics_by_replicate.csv"
  ),
  row.names = FALSE
)

write.csv(
  scores_all_nodes,
  file.path(
    OUT_DIR,
    "competing_methods_all_node_probabilities.csv"
  ),
  row.names = FALSE
)

write.csv(
  geometry_summary,
  file.path(
    OUT_DIR,
    "simulation_geometry_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  modularity_summary,
  file.path(
    OUT_DIR,
    "simulation_modularity_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  failure_log,
  file.path(
    OUT_DIR,
    "competing_methods_failure_log.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 13. MONTE CARLO SUMMARIES
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
  Overall_Coverage_80 = "higher_interpret_with_interval_score",
  Hidden_Coverage_80 = "higher_interpret_with_interval_score",
  Overall_Interval_Score_80 = "lower",
  Hidden_Interval_Score_80 = "lower",
  Overall_Coverage_95 = "higher_interpret_with_interval_score",
  Hidden_Coverage_95 = "higher_interpret_with_interval_score",
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
      "Method ", method_name,
      " does not have exactly 100 replicate rows."
    )
  }
  
  for (metric_name in main_metric_names) {
    statistics <- summary_stats(method_data[[metric_name]])
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

metric_summary_long <- do.call(rbind, summary_rows)
row.names(metric_summary_long) <- NULL

write.csv(
  metric_summary_long,
  file.path(
    OUT_DIR,
    "competing_methods_metric_summary_long.csv"
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
    "competing_methods_FINAL_MAIN_TABLE_means.csv"
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
    return(formatC(mean_value, format = "f", digits = digits))
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
        if (!is.finite(x)) "" else formatC(x, format = "f", digits = 3L)
      },
      character(1)
    )
  } else {
    digits_now <- 3L
    formatted <- mapply(
      format_mean_sd,
      metric_rows$mean,
      metric_rows$sd,
      MoreArgs = list(digits = digits_now),
      USE.NAMES = FALSE
    )
  }
  
  paper_table_mean_sd[[paper_column_labels[metric_name]]] <- formatted
}


## ============================================================================
## ADD MONTE CARLO SD TO ALL COVERAGE COLUMNS
## ============================================================================

coverage_metric_to_column <- c(
  Overall_Coverage_80 = "Overall Coverage 80%",
  Hidden_Coverage_80  = "Hidden Coverage 80%",
  Overall_Coverage_95 = "Overall Coverage 95%",
  Hidden_Coverage_95  = "Hidden Coverage 95%"
)

format_coverage_mean_sd <- function(values, digits = 3L) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  
  ## Ada, Ada + Net and GNN have no uncertainty results.
  if (length(values) == 0L) {
    return("")
  }
  
  paste0(
    formatC(
      mean(values),
      format = "f",
      digits = digits
    ),
    " (",
    formatC(
      sd(values),
      format = "f",
      digits = digits
    ),
    ")"
  )
}

for (metric_name in names(coverage_metric_to_column)) {
  
  output_column <- unname(
    coverage_metric_to_column[metric_name]
  )
  
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

## Rewrite the paper table.
write.csv(
  paper_table_mean_sd,
  file.path(
    OUT_DIR,
    "competing_methods_FINAL_PAPER_TABLE_mean_SD.csv"
  ),
  row.names = FALSE
)

cat(
  "\n================ FINAL PAPER TABLE WITH COVERAGE MEAN (SD) ================\n"
)

print(
  paper_table_mean_sd,
  row.names = FALSE
)

write.csv(
  paper_table_mean_sd,
  file.path(
    OUT_DIR,
    "competing_methods_FINAL_PAPER_TABLE_mean_SD.csv"
  ),
  row.names = FALSE,
  na = ""
)

## A compact uncertainty-only table is also saved for convenient inspection.
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
    unname(paper_column_labels[uncertainty_metric_names])
  ),
  drop = FALSE
]

write.csv(
  uncertainty_table_mean_sd,
  file.path(
    OUT_DIR,
    "competing_methods_uncertainty_metrics_80_95.csv"
  ),
  row.names = FALSE,
  na = ""
)


## ============================================================================
## 14. SAVE COMPLETE OUTPUT BUNDLE
## ============================================================================

total_runtime_sec <- as.numeric(
  difftime(
    simulation_end,
    simulation_start,
    units = "secs"
  )
)

run_config <- list(
  script_version = SCRIPT_VERSION,
  created_at = timestamp_now(),
  scenario = SCENARIO_NAME,
  settings_signature = SETTINGS_SIGNATURE,
  output_directory = OUT_DIR,
  n_rep = N_REP,
  n_nodes = N_NODES,
  true_labels = TRUE_T,
  observed_labels = OBSERVED_Y,
  signal_type = SIGNAL_TYPE,
  true_zero_idx = TRUE_ZERO_IDX,
  hidden_positive_idx = HIDDEN_POS_IDX,
  observed_positive_idx = OBSERVED_POS_IDX,
  hidden_covariate_only_idx = HIDDEN_COV_ONLY_IDX,
  hidden_network_only_idx = HIDDEN_NET_ONLY_IDX,
  hidden_both_idx = HIDDEN_BOTH_IDX,
  observed_covariate_only_idx = OBSERVED_COV_ONLY_IDX,
  observed_network_only_idx = OBSERVED_NET_ONLY_IDX,
  observed_both_idx = OBSERVED_BOTH_IDX,
  covariate_signal_idx = COVARIATE_SIGNAL_IDX,
  network_signal_idx = NETWORK_SIGNAL_IDX,
  network_block = NETWORK_BLOCK,
  covariate_dgp = list(
    n_covariates = N_COVARIATES,
    background_mean = COVARIATE_BACKGROUND_MEAN,
    signal_mean = COVARIATE_SIGNAL_MEAN,
    common_sd = COMMON_SD,
    allocation = "covariate-only plus both are signal"
  ),
  network_dgp = list(
    type = "two-block Bernoulli SBM",
    block_definition = "network-only plus both form block 1",
    p_00 = P_00,
    p_01 = P_01,
    p_11 = P_11,
    triad_closure_prob = TRIAD_CLOSURE_PROB,
    expected_degree_network_background =
      EXPECTED_DEGREE_NETWORK_BACKGROUND,
    expected_degree_network_signal =
      EXPECTED_DEGREE_NETWORK_SIGNAL
  ),
  methods = method_registry,
  metric_definitions = list(
    Overall_AUC = paste0(
      "AUC on all 400 observations against true T with 45 true positives."
    ),
    Hidden_AUC = paste0(
      "AUC among the 370 Y=0 observations, containing 15 hidden positives ",
      "and 355 true zeros."
    ),
    Overall_LogLoss = "LogLoss on all 400 observations against true T.",
    Hidden_LogLoss = paste0(
      "Mean -log(p_i) over the 15 hidden positives, all of which have T_i=1."
    ),
    Overall_Brier = "Brier score on all 400 observations against true T.",
    Hidden_Brier = paste0(
      "Mean (p_i-1)^2 over the 15 hidden positives."
    ),
    Overall_Recall_at_45 = paste0(
      "Fraction of all 45 true positives whose descending midrank among ",
      "all 400 observations is at most 45."
    ),
    Hidden_Recall_at_15 = paste0(
      "Fraction of the 15 hidden positives whose descending midrank among ",
      "the 370 Y=0 observations is at most 15."
    ),
    Coverage = paste0(
      "For each node, the central empirical binary predictive interval is ",
      "formed from predictive T draws using quantile type 1 and compared ",
      "with known TRUE_T. Overall averages over 400 nodes; Hidden averages ",
      "over the 15 hidden positives."
    ),
    Interval_Score = paste0(
      "Equal-tail interval score evaluated against known TRUE_T, with ",
      "alpha=0.20 at 80% and alpha=0.05 at 95%."
    )
  ),
  ranking_rule = "rank(-probability, ties.method='average')",
  predictive_uncertainty = list(
    levels = PREDICTIVE_INTERVAL_LEVELS,
    draws = N_PREDICTIVE_DRAWS,
    empirical_quantile_type = 1L,
    WLR = paste0(
      "Asymptotic Gaussian uncertainty for each linear predictor using ",
      "x_i' vcov(beta_hat) x_i, followed by Bernoulli predictive draws."
    ),
    INLA = paste0(
      "Posterior linear-predictor marginal draws from INLA, followed by ",
      "Bernoulli predictive draws."
    ),
    Ada = "NA: no model-based uncertainty in the implemented method.",
    Ada_Net = "NA: no model-based uncertainty in the implemented method.",
    GNN = "NA: no model-based uncertainty in the implemented method.",
    poisson_binomial = "not used",
    H_coverage = "not computed"
  ),
  seeds = seed_table,
  resume_replicates = RESUME_REPLICATES,
  force_fresh_run = FORCE_FRESH_RUN,
  n_loaded_complete = n_loaded_complete,
  n_computed_now = n_computed_now,
  total_runtime_sec = total_runtime_sec
)

saveRDS(
  run_config,
  file.path(
    OUT_DIR,
    "simulation_run_config.rds"
  )
)

save(
  replication_results,
  metrics_by_replicate,
  scores_all_nodes,
  geometry_summary,
  modularity_summary,
  failure_log,
  metric_summary_long,
  main_table_means,
  paper_table_mean_sd,
  uncertainty_table_mean_sd,
  method_registry,
  seed_table,
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
  ZERO_MEAN,
  POSITIVE_MEAN,
  COMMON_SD,
  P_00,
  P_01,
  P_11,
  TRIAD_CLOSURE_PROB,
  run_config,
  SETTINGS_SIGNATURE,
  file = file.path(
    OUT_DIR,
    "SIMULATION_100_COMPETING_METHODS_LINEAR_GAUSSIAN_SBM_SCAR_MIXED_SIGNAL_COMPLETE_RESULTS.RData"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    OUT_DIR,
    "sessionInfo_competing_methods_linear_gaussian_SBM_SCAR_mixed_signal.txt"
  )
)


## ============================================================================
## 15. FINAL CONSOLE OUTPUT
## ============================================================================

cat("\n================ RUN ACCOUNTING ================\n")
cat("Requested replications :", N_REP, "\n")
cat("Computed in this run   :", n_computed_now, "\n")
cat("Loaded from checkpoint :", n_loaded_complete, "\n")
cat("Total runtime          :", format_duration(total_runtime_sec), "\n")

cat("\n================ FINAL MONTE CARLO MEANS ================\n")
print(
  main_table_means,
  row.names = FALSE
)

cat("\n================ FINAL PAPER TABLE: MEAN (SD) ================\n")
print(
  paper_table_mean_sd,
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
