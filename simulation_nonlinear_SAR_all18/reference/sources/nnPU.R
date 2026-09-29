## ============================================================================
## SIMULATION STUDY: nnPU ON THE SAME 100 MIXED-SIGNAL SAR DATA SETS
## ----------------------------------------------------------------------------
## This script fits the non-negative Positive-Unlabeled method of
## Kiryo et al. (2017) in two versions:
##
##   1) nnPU       : the ten standardized covariates only;
##   2) nnPU + Net : the same covariates plus the label-blind leading
##                  eigenvectors of the modularity matrix.
##
## DATA AND EVALUATION
## -------------------
## The latent 400-node mixed-signal geometry is the same as in the nonlinear
## SCAR study. The only substantive data change is the label mechanism:
##
##   - 45 true positives are retained in every replication;
##   - exactly 15 are hidden under the covariate-driven SAR mechanism with
##     lambda = 2 used by the competing-method SAR simulation;
##   - the remaining 30 positives are observed;
##   - the SAR oracle score is used only to generate missing labels and is never
##     supplied to nnPU.
##
## Production requires and loads the exact saved SAR geometry from
##
##   SIMULATION_100_COMPETING_METHODS_MIXED_SIGNAL_SAR_LAMBDA_2/
##     generated_geometries
##
## after strict validation of labels, seeds, DGP settings, SAR settings and
## modularity features. Therefore every nnPU fit is paired with the same
## X, A, Z, observed-positive set and hidden-positive set as the competing
## methods.
##
## nnPU IMPLEMENTATION
## -------------------
## The classifier, risk, optimizer, training length and evaluation are unchanged
## from the SCAR nnPU script:
##
##   g(x) = W2' ReLU(W1' x + b1) + b2,
##
## with 100 hidden units, sigmoid surrogate loss, Chainer-style Adam, L2 weight
## decay 0.005, beta = 0, gamma = 1 and 5,000 full-batch updates.
##
## In every SAR replication:
##
##   P = the 30 replication-specific observed positive rows;
##   U = all 400 nodes, representing the marginal p(x).
##
## Hence every update uses 430 training rows. The class prior remains
##
##   pi_p = 45 / 400 = 0.1125.
##
## AUC and recall use raw decision scores. LogLoss and Brier use plogis(score)
## without post-hoc calibration. Predictive coverage and interval-score columns
## are retained as NA because the official nnPU implementation is a
## point-estimation classifier.
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

SCRIPT_VERSION <- "simulation_nnpu_mixed_signal_SAR_lambda2_same_data_v1_0"
MASTER_SEED <- 20260718L
N_REP <- 100L
N_NODES <- 400L

## --------------------------------------------------------------------------
## 1A. Fixed latent truth and fixed mixed-signal geometry
## --------------------------------------------------------------------------

TRUE_ZERO_IDX <- 1:355
TRUE_POS_IDX <- 356:400

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

N_NNPU_MARGINAL_U <- N_NODES
N_NNPU_TRAIN_ROWS <- N_OBSERVED_POS + N_NNPU_MARGINAL_U

TRUE_T <- integer(N_NODES)
TRUE_T[TRUE_POS_IDX] <- 1L

POSITIVE_SIGNAL_GROUP <- rep(
  "true zero",
  N_NODES
)
POSITIVE_SIGNAL_GROUP[COV_ONLY_IDX] <- "covariates only"
POSITIVE_SIGNAL_GROUP[NET_ONLY_IDX] <- "network only"
POSITIVE_SIGNAL_GROUP[BOTH_IDX] <- "both"

## --------------------------------------------------------------------------
## 1B. SAR label mechanism
## --------------------------------------------------------------------------

LAMBDA_SAR <- 2.00
SAR_Z_CLIP <- 3.00

SCENARIO_NAME <- paste0(
  "Fixed mixed-signal covariate-driven SAR design, lambda = ",
  format(LAMBDA_SAR, trim = TRUE)
)

## --------------------------------------------------------------------------
## 1C. nnPU settings: unchanged from the SCAR script
## --------------------------------------------------------------------------

NNPU_CLASS_PRIOR <- 45 / 400
NNPU_HIDDEN_UNITS <- 100L
NNPU_EPOCHS <- 5000L
NNPU_STEPSIZE <- 1e-3
NNPU_WEIGHT_DECAY <- 0.005
NNPU_BETA <- 0
NNPU_GAMMA <- 1
NNPU_ADAM_BETA1 <- 0.9
NNPU_ADAM_BETA2 <- 0.999
NNPU_ADAM_EPS <- 1e-8
NNPU_PROBABILITY_MAP <- "plogis(raw decision score); no post-hoc calibration"

PREDICTIVE_INTERVAL_LEVELS <- c(0.80, 0.95)
N_PREDICTIVE_DRAWS <- 2000L

HIDDEN_TOP_K <- 15L
OVERALL_RECALL_K <- 45L
PROBABILITY_EPS <- 1e-8

## --------------------------------------------------------------------------
## 1D. Covariate DGP settings
## --------------------------------------------------------------------------

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

## --------------------------------------------------------------------------
## 1E. Network DGP settings
## --------------------------------------------------------------------------

ALPHA_BACKGROUND <- 2.0
ALPHA_FRAUD <- 2.0
ALPHA_LEGAL <- 2.0

D_BACKGROUND <- 2
D_FRAUD <- 2
D_LEGAL <- 2

TRIAD_CLOSURE_PROB <- 0.05
TRIAD_MAX_NEW_EDGES <- Inf

## --------------------------------------------------------------------------
## 1F. Modularity-feature settings
## --------------------------------------------------------------------------

MOD_EIG_TOL <- 1e-8
MOD_MAX_Q <- Inf
MOD_K_MAX <- 20L

## --------------------------------------------------------------------------
## 1G. Output, restart and exact SAR geometry reuse
## --------------------------------------------------------------------------

OUT_DIR <- file.path(
  getwd(),
  "SIMULATION_100_NNPU_MIXED_SIGNAL_SAR_LAMBDA_2"
)
CHECKPOINT_DIR <- file.path(OUT_DIR, "replicate_checkpoints")
GEOMETRY_DIR <- file.path(OUT_DIR, "generated_geometries")

for (path_now in c(OUT_DIR, CHECKPOINT_DIR, GEOMETRY_DIR)) {
  dir.create(path_now, showWarnings = FALSE, recursive = TRUE)
}

RESUME_REPLICATES <- TRUE
FORCE_FRESH_RUN <- FALSE
FAIL_IF_ANY_METHOD_FAILS <- TRUE
SAVE_GENERATED_GEOMETRIES <- TRUE
SAVE_TRAINING_TRACES <- TRUE

USE_SAVED_COMPETING_GEOMETRIES <- TRUE
REQUIRE_SAVED_COMPETING_GEOMETRIES <- TRUE

COMPETING_OUT_DIR <- file.path(
  getwd(),
  "SIMULATION_100_COMPETING_METHODS_MIXED_SIGNAL_SAR_LAMBDA_2"
)
COMPETING_GEOMETRY_DIR <- file.path(
  COMPETING_OUT_DIR,
  "generated_geometries"
)

MERGE_WITH_EXISTING_COMPETITOR_RESULTS <- TRUE

## --------------------------------------------------------------------------
## 1H. Deterministic seeds
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
  nnpu_seed = as.integer(
    MASTER_SEED + 1200000L + 1000L * seq_len(N_REP) + 149L
  ),
  nnpu_net_seed = as.integer(
    MASTER_SEED + 1200000L + 1000L * seq_len(N_REP) + 149L
  ),
  stringsAsFactors = FALSE
)

seed_columns <- setdiff(names(seed_table), "replicate")
if (any(as.matrix(seed_table[, seed_columns, drop = FALSE]) >
        .Machine$integer.max)) {
  stop("At least one deterministic seed exceeds the R integer range.")
}

write.csv(
  seed_table,
  file.path(OUT_DIR, "nnpu_SAR_simulation_seed_table.csv"),
  row.names = FALSE
)

## ============================================================================
## 2. DESIGN CHECKS
## ============================================================================

if (N_REP != 100L) stop("The requested simulation study requires N_REP = 100.")
if (N_NODES != 400L) stop("The requested simulation study requires 400 nodes.")
if (!identical(TRUE_ZERO_IDX, 1:355)) stop("TRUE_ZERO_IDX must be 1:355.")
if (!identical(TRUE_POS_IDX, 356:400)) stop("TRUE_POS_IDX must be 356:400.")

if (
  length(COV_ONLY_IDX) != 15L ||
  length(NET_ONLY_IDX) != 15L ||
  length(BOTH_IDX) != 15L ||
  length(unique(c(COV_ONLY_IDX, NET_ONLY_IDX, BOTH_IDX))) != N_TRUE_POS ||
  !identical(
    sort(c(COV_ONLY_IDX, NET_ONLY_IDX, BOTH_IDX)),
    TRUE_POS_IDX
  )
) {
  stop("The three positive signal groups must partition TRUE_POS_IDX.")
}

if (
  N_TRUE_POS != 45L ||
  N_HIDDEN != 15L ||
  N_OBSERVED_POS != 30L ||
  N_TRUE_ZERO != 355L ||
  N_UNLABELED != 370L ||
  N_NNPU_MARGINAL_U != 400L ||
  N_NNPU_TRAIN_ROWS != 430L
) {
  stop("The fixed SAR PU design or the 30-P/400-U nnPU construction changed.")
}

if (
  length(COVARIATE_SIGNAL_IDX) != 30L ||
  length(NETWORK_SIGNAL_IDX) != 30L ||
  length(intersect(COVARIATE_SIGNAL_IDX, NETWORK_SIGNAL_IDX)) != 15L
) {
  stop("The mixed-signal geometry is inconsistent.")
}

if (!identical(which(TRUE_T == 1L), TRUE_POS_IDX)) {
  stop("TRUE_T is inconsistent with TRUE_POS_IDX.")
}

if (
  !isTRUE(all.equal(LAMBDA_SAR, 2.00)) ||
  !is.finite(SAR_Z_CLIP) ||
  SAR_Z_CLIP <= 0
) {
  stop("The requested SAR mechanism requires lambda=2 and positive clipping.")
}

if (!isTRUE(all.equal(NNPU_CLASS_PRIOR, 0.1125))) {
  stop("NNPU_CLASS_PRIOR must equal 45/400 = 0.1125.")
}

if (
  NNPU_HIDDEN_UNITS != 100L ||
  NNPU_EPOCHS != 5000L ||
  NNPU_STEPSIZE <= 0 ||
  NNPU_WEIGHT_DECAY < 0 ||
  NNPU_BETA < 0 ||
  NNPU_GAMMA < 0 ||
  NNPU_GAMMA > 1
) {
  stop("Invalid nnPU hyperparameters.")
}

if (
  HIDDEN_TOP_K != 15L ||
  OVERALL_RECALL_K != 45L ||
  !isTRUE(all.equal(PREDICTIVE_INTERVAL_LEVELS, c(0.80, 0.95)))
) {
  stop("The requested evaluation settings have changed.")
}

if (
  isTRUE(REQUIRE_SAVED_COMPETING_GEOMETRIES) &&
  !isTRUE(USE_SAVED_COMPETING_GEOMETRIES)
) {
  stop(
    "REQUIRE_SAVED_COMPETING_GEOMETRIES=TRUE requires ",
    "USE_SAVED_COMPETING_GEOMETRIES=TRUE."
  )
}

if (isTRUE(FORCE_FRESH_RUN)) {
  unlink(CHECKPOINT_DIR, recursive = TRUE, force = TRUE)
  dir.create(CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)
}

cat("\n================ nnPU SAR SIMULATION DESIGN ================\n")
cat("Replications             :", N_REP, "\n")
cat("Observations             :", N_NODES, "\n")
cat("True positives           :", N_TRUE_POS, "\n")
cat("Observed positives       :", N_OBSERVED_POS, "per replication\n")
cat("Hidden positives         :", N_HIDDEN, "per replication\n")
cat("Y=0 evaluation pool      :", N_UNLABELED, "per replication\n")
cat("SAR lambda               :", LAMBDA_SAR, "\n")
cat("nnPU training P / U      :", N_OBSERVED_POS, "/", N_NNPU_MARGINAL_U, "\n")
cat("nnPU training rows       :", N_NNPU_TRAIN_ROWS, "\n")
cat("Class prior pi_p         :", NNPU_CLASS_PRIOR, "\n")
cat("Model                    : d-100-1 ReLU\n")
cat("Full-batch updates / lr  :", NNPU_EPOCHS, "/", NNPU_STEPSIZE, "\n")
cat("beta / gamma             :", NNPU_BETA, "/", NNPU_GAMMA, "\n")
cat("Weight decay             :", NNPU_WEIGHT_DECAY, "\n")

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
    hidden_idx,
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
    hidden_idx
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


## ============================================================================
## 7. OFFICIAL-STYLE FULL-BATCH nnPU 3LP IN BASE R
## ============================================================================

## Chainer's default Linear initialization is a LeCun normal initialization.
initialize_nnpu_3lp <- function(input_dimension, hidden_units, seed) {
  input_dimension <- as.integer(input_dimension)
  hidden_units <- as.integer(hidden_units)
  if (input_dimension < 1L || hidden_units < 1L) {
    stop("input_dimension and hidden_units must be positive.")
  }

  set.seed(as.integer(seed))

  list(
    W1 = matrix(
      rnorm(
        input_dimension * hidden_units,
        mean = 0,
        sd = sqrt(1 / input_dimension)
      ),
      nrow = input_dimension,
      ncol = hidden_units
    ),
    b1 = rep(0, hidden_units),
    W2 = matrix(
      rnorm(
        hidden_units,
        mean = 0,
        sd = sqrt(1 / hidden_units)
      ),
      nrow = hidden_units,
      ncol = 1L
    ),
    b2 = 0
  )
}

zero_like_parameter_list <- function(parameters) {
  lapply(parameters, function(x) x * 0)
}

nnpu_forward_3lp <- function(X, parameters) {
  X <- as.matrix(X)
  preactivation <- sweep(
    X %*% parameters$W1,
    2L,
    parameters$b1,
    "+"
  )
  hidden <- pmax(preactivation, 0)
  score <- as.vector(hidden %*% parameters$W2 + parameters$b2)

  list(
    preactivation = preactivation,
    hidden = hidden,
    score = score
  )
}

nnpu_risk_and_score_gradient <- function(
    score,
    pu_target,
    class_prior = NNPU_CLASS_PRIOR,
    beta = NNPU_BETA,
    gamma = NNPU_GAMMA
) {
  score <- as.numeric(score)
  pu_target <- as.integer(pu_target)

  positive <- pu_target == 1L
  unlabeled <- pu_target == -1L
  n_positive <- sum(positive)
  n_unlabeled <- sum(unlabeled)

  if (n_positive < 1L || n_unlabeled < 1L) {
    stop("Every nnPU full batch must contain at least one P and one U case.")
  }
  if (!is.finite(class_prior) || class_prior <= 0 || class_prior >= 1) {
    stop("class_prior must lie in (0,1).")
  }

  probability <- plogis(score)
  sigmoid_derivative <- probability * (1 - probability)

  ## Authors' sigmoid loss:
  ##   loss(score,+1) = sigmoid(-score) = 1 - sigmoid(score)
  ##   loss(score,-1) = sigmoid( score)
  positive_risk <- class_prior * mean(1 - probability[positive])
  negative_risk <- mean(probability[unlabeled]) -
    class_prior * mean(probability[positive])

  objective <- positive_risk + negative_risk
  correction_active <- isTRUE(negative_risk < -beta)

  if (correction_active) {
    ## Exact pu_loss.py correction gradient: differentiate -gamma * R_N.
    differentiated_objective <- -gamma * negative_risk
    score_gradient <- -gamma * (
      as.numeric(unlabeled) / n_unlabeled -
        class_prior * as.numeric(positive) / n_positive
    ) * sigmoid_derivative

    ## This is what the official forward method returns for logging.
    reported_objective <- positive_risk - beta
  } else {
    differentiated_objective <- objective

    ## Gradient of pi*R_P+ + R_U- - pi*R_P-.
    score_gradient <- (
      as.numeric(unlabeled) / n_unlabeled -
        2 * class_prior * as.numeric(positive) / n_positive
    ) * sigmoid_derivative

    reported_objective <- objective
  }

  list(
    positive_risk = as.numeric(positive_risk),
    negative_risk = as.numeric(negative_risk),
    unbiased_objective = as.numeric(objective),
    nonnegative_objective = as.numeric(
      positive_risk + max(0, negative_risk)
    ),
    reported_objective = as.numeric(reported_objective),
    differentiated_objective = as.numeric(differentiated_objective),
    correction_active = correction_active,
    score_gradient = as.numeric(score_gradient)
  )
}

nnpu_parameter_gradients <- function(X, forward, score_gradient, parameters) {
  X <- as.matrix(X)
  score_gradient <- as.numeric(score_gradient)

  gradient_W2 <- crossprod(forward$hidden, score_gradient)
  gradient_b2 <- sum(score_gradient)

  gradient_hidden <- tcrossprod(
    score_gradient,
    as.vector(parameters$W2)
  )
  gradient_preactivation <- gradient_hidden *
    (forward$preactivation > 0)

  gradient_W1 <- crossprod(X, gradient_preactivation)
  gradient_b1 <- colSums(gradient_preactivation)

  list(
    W1 = gradient_W1,
    b1 = gradient_b1,
    W2 = matrix(gradient_W2, ncol = 1L),
    b2 = as.numeric(gradient_b2)
  )
}

add_weight_decay_to_gradients <- function(
    gradients,
    parameters,
    weight_decay
) {
  if (!is.finite(weight_decay) || weight_decay < 0) {
    stop("weight_decay must be non-negative.")
  }

  for (parameter_name in names(parameters)) {
    gradients[[parameter_name]] <- gradients[[parameter_name]] +
      weight_decay * parameters[[parameter_name]]
  }
  gradients
}

## Chainer-style Adam: bias correction is absorbed into the time-varying
## learning rate, while epsilon is added to sqrt(v) before the update.
adam_parameter_step <- function(
    parameter,
    gradient,
    first_moment,
    second_moment,
    iteration,
    stepsize,
    beta1,
    beta2,
    epsilon
) {
  first_moment <- beta1 * first_moment +
    (1 - beta1) * gradient
  second_moment <- beta2 * second_moment +
    (1 - beta2) * (gradient * gradient)

  learning_rate_now <- stepsize *
    sqrt(1 - beta2^iteration) /
    (1 - beta1^iteration)

  parameter <- parameter -
    learning_rate_now * first_moment /
    (sqrt(second_moment) + epsilon)

  list(
    parameter = parameter,
    first_moment = first_moment,
    second_moment = second_moment
  )
}

fit_nnpu_3lp <- function(
    X,
    observed_y,
    seed,
    class_prior = NNPU_CLASS_PRIOR,
    hidden_units = NNPU_HIDDEN_UNITS,
    epochs = NNPU_EPOCHS,
    stepsize = NNPU_STEPSIZE,
    weight_decay = NNPU_WEIGHT_DECAY,
    beta = NNPU_BETA,
    gamma = NNPU_GAMMA,
    adam_beta1 = NNPU_ADAM_BETA1,
    adam_beta2 = NNPU_ADAM_BETA2,
    adam_eps = NNPU_ADAM_EPS
) {
  X <- as.matrix(X)
  observed_y <- as.integer(observed_y)
  n <- nrow(X)
  d <- ncol(X)

  if (n != length(observed_y) || n < 2L || d < 1L) {
    stop("X and observed_y have incompatible dimensions.")
  }
  if (any(!is.finite(X)) || !all(observed_y %in% c(0L, 1L))) {
    stop("nnPU received non-finite features or invalid observed labels.")
  }

  ## nnPU two-sample construction:
  ##   P = the 30 observed positive rows;
  ##   U = all 400 nodes, representing the marginal p(x).
  ##
  ## The 30 observed positives therefore occur once as P and once inside U.
  positive_idx <- which(observed_y == 1L)

  X_train <- rbind(
    X[positive_idx, , drop = FALSE],
    X
  )

  ## Official PULoss labels: +1 = labeled positive, -1 = unlabeled.
  pu_target <- c(
    rep.int(1L, length(positive_idx)),
    rep.int(-1L, n)
  )
  storage.mode(pu_target) <- "integer"

  if (
    sum(pu_target == 1L) != N_OBSERVED_POS ||
    sum(pu_target == -1L) != N_NNPU_MARGINAL_U ||
    nrow(X_train) != N_NNPU_TRAIN_ROWS
  ) {
    stop(
      "The nnPU training construction must contain 30 P rows, ",
      "400 marginal-U rows and 430 rows in total."
    )
  }

  parameters <- initialize_nnpu_3lp(
    input_dimension = d,
    hidden_units = hidden_units,
    seed = seed
  )
  first_moment <- zero_like_parameter_list(parameters)
  second_moment <- zero_like_parameter_list(parameters)

  trace <- data.frame(
    epoch = seq_len(as.integer(epochs)),
    positive_risk = NA_real_,
    negative_risk = NA_real_,
    unbiased_objective = NA_real_,
    nonnegative_objective = NA_real_,
    reported_objective = NA_real_,
    differentiated_objective = NA_real_,
    correction_active = NA_integer_,
    max_abs_score = NA_real_,
    parameter_l2_norm = NA_real_,
    stringsAsFactors = FALSE
  )

  fit_start <- Sys.time()

  for (epoch in seq_len(as.integer(epochs))) {
    forward <- nnpu_forward_3lp(X_train, parameters)
    risk <- nnpu_risk_and_score_gradient(
      score = forward$score,
      pu_target = pu_target,
      class_prior = class_prior,
      beta = beta,
      gamma = gamma
    )

    gradients <- nnpu_parameter_gradients(
      X = X_train,
      forward = forward,
      score_gradient = risk$score_gradient,
      parameters = parameters
    )
    gradients <- add_weight_decay_to_gradients(
      gradients = gradients,
      parameters = parameters,
      weight_decay = weight_decay
    )

    for (parameter_name in names(parameters)) {
      update <- adam_parameter_step(
        parameter = parameters[[parameter_name]],
        gradient = gradients[[parameter_name]],
        first_moment = first_moment[[parameter_name]],
        second_moment = second_moment[[parameter_name]],
        iteration = epoch,
        stepsize = stepsize,
        beta1 = adam_beta1,
        beta2 = adam_beta2,
        epsilon = adam_eps
      )

      parameters[[parameter_name]] <- update$parameter
      first_moment[[parameter_name]] <- update$first_moment
      second_moment[[parameter_name]] <- update$second_moment
    }

    trace$positive_risk[epoch] <- risk$positive_risk
    trace$negative_risk[epoch] <- risk$negative_risk
    trace$unbiased_objective[epoch] <- risk$unbiased_objective
    trace$nonnegative_objective[epoch] <- risk$nonnegative_objective
    trace$reported_objective[epoch] <- risk$reported_objective
    trace$differentiated_objective[epoch] <-
      risk$differentiated_objective
    trace$correction_active[epoch] <- as.integer(
      risk$correction_active
    )
    trace$max_abs_score[epoch] <- max(abs(forward$score))
    trace$parameter_l2_norm[epoch] <- sqrt(sum(
      unlist(parameters, use.names = FALSE)^2
    ))

    if (epoch %% 10L == 0L) {
      parameter_vector <- unlist(parameters, use.names = FALSE)
      if (any(!is.finite(parameter_vector))) {
        stop("nnPU parameters became non-finite at epoch ", epoch, ".")
      }
    }
  }

  ## Evaluate the final nnPU training risk on the 430-row P/U construction.
  final_train_forward <- nnpu_forward_3lp(
    X_train,
    parameters
  )
  final_risk <- nnpu_risk_and_score_gradient(
    score = final_train_forward$score,
    pu_target = pu_target,
    class_prior = class_prior,
    beta = beta,
    gamma = gamma
  )

  ## Produce one score and probability for each of the 400 unique nodes.
  final_prediction_forward <- nnpu_forward_3lp(
    X,
    parameters
  )
  probabilities <- plogis(final_prediction_forward$score)

  if (any(!is.finite(final_prediction_forward$score)) ||
      any(!is.finite(probabilities))) {
    stop("nnPU returned non-finite final scores or probabilities.")
  }

  runtime_sec <- as.numeric(difftime(
    Sys.time(),
    fit_start,
    units = "secs"
  ))

  diagnostics <- data.frame(
    seed = as.integer(seed),
    n = n,
    n_unique_nodes = n,
    n_training_rows = nrow(X_train),
    input_dimension = d,
    n_positive = sum(pu_target == 1L),
    n_unlabeled = sum(pu_target == -1L),
    class_prior = class_prior,
    hidden_units = as.integer(hidden_units),
    epochs = as.integer(epochs),
    stepsize = stepsize,
    weight_decay = weight_decay,
    beta = beta,
    gamma = gamma,
    final_positive_risk = final_risk$positive_risk,
    final_negative_risk = final_risk$negative_risk,
    final_unbiased_objective = final_risk$unbiased_objective,
    final_nonnegative_objective = final_risk$nonnegative_objective,
    final_reported_objective = final_risk$reported_objective,
    final_correction_active = as.integer(final_risk$correction_active),
    correction_steps = sum(trace$correction_active),
    correction_fraction = mean(trace$correction_active),
    score_min = min(final_prediction_forward$score),
    score_mean = mean(final_prediction_forward$score),
    score_sd = stats::sd(final_prediction_forward$score),
    score_max = max(final_prediction_forward$score),
    probability_min = min(probabilities),
    probability_mean = mean(probabilities),
    probability_max = max(probabilities),
    parameter_l2_norm = sqrt(sum(
      unlist(parameters, use.names = FALSE)^2
    )),
    internal_fit_runtime_sec = runtime_sec,
    stringsAsFactors = FALSE
  )

  list(
    scores = as.numeric(final_prediction_forward$score),
    probabilities = as.numeric(probabilities),
    parameters = parameters,
    trace = trace,
    diagnostics = diagnostics
  )
}

method_registry <- data.frame(
  method_id = c("NNPU", "NNPU_NET"),
  method = c("nnPU", "nnPU + Net"),
  stringsAsFactors = FALSE
)

run_nnpu_method <- function(
    method_id,
    X_cov,
    X_cov_net,
    y,
    replicate_seed_row
) {
  if (identical(method_id, "NNPU")) {
    X_fit <- X_cov
    fit_seed <- replicate_seed_row$nnpu_seed
  } else if (identical(method_id, "NNPU_NET")) {
    X_fit <- X_cov_net
    fit_seed <- replicate_seed_row$nnpu_net_seed
  } else {
    stop("Unknown nnPU method_id: ", method_id)
  }

  fit <- fit_nnpu_3lp(
    X = X_fit,
    observed_y = y,
    seed = fit_seed,
    class_prior = NNPU_CLASS_PRIOR,
    hidden_units = NNPU_HIDDEN_UNITS,
    epochs = NNPU_EPOCHS,
    stepsize = NNPU_STEPSIZE,
    weight_decay = NNPU_WEIGHT_DECAY,
    beta = NNPU_BETA,
    gamma = NNPU_GAMMA
  )

  list(
    probabilities = fit$probabilities,
    scores = fit$scores,
    predictive_T_draws = NULL,
    uncertainty_type = paste0(
      "not available: official nnPU is a point-estimation classifier; ",
      "probabilities are plogis(raw score)"
    ),
    training_trace = fit$trace,
    training_diagnostics = fit$diagnostics
  )
}


## ============================================================================
## 8. COMMON EVALUATION
## ============================================================================

evaluate_method_scores <- function(
    method_output,
    method_id,
    method_label,
    replicate_number,
    observed_y,
    hidden_positive_idx,
    signal_type
) {
  if (!is.list(method_output) || is.null(method_output$probabilities)) {
    stop(method_label, " did not return the required method-output list.")
  }

  observed_y <- as.integer(observed_y)
  hidden_positive_idx <- sort(unique(as.integer(hidden_positive_idx)))
  signal_type <- as.character(signal_type)

  if (
    length(observed_y) != N_NODES ||
    !all(observed_y %in% c(0L, 1L)) ||
    sum(observed_y) != N_OBSERVED_POS
  ) {
    stop("The replication-specific SAR observed-label vector is invalid.")
  }

  if (
    length(hidden_positive_idx) != N_HIDDEN ||
    any(!hidden_positive_idx %in% TRUE_POS_IDX) ||
    any(observed_y[hidden_positive_idx] != 0L)
  ) {
    stop("The replication-specific SAR hidden-positive indices are invalid.")
  }

  if (length(signal_type) != N_NODES || any(is.na(signal_type))) {
    stop("The replication-specific signal_type vector is invalid.")
  }

  probabilities <- validate_probability_vector(
    probabilities = method_output$probabilities,
    n_expected = N_NODES,
    method_label = method_label
  )

  decision_scores <- as.numeric(method_output$scores)
  if (
    length(decision_scores) != N_NODES ||
    any(!is.finite(decision_scores))
  ) {
    stop(method_label, " returned invalid raw decision scores.")
  }

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
    -decision_scores,
    ties.method = "average"
  )

  unlabeled_probabilities <- probabilities[unlabeled_idx]
  unlabeled_decision_scores <- decision_scores[unlabeled_idx]

  unlabeled_rank <- rank(
    -unlabeled_decision_scores,
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
      hidden_idx = hidden_positive_idx,
      level = 0.80
    )

    predictive_95 <- predictive_label_metrics(
      predictive_T_draws = predictive_T_draws,
      truth = TRUE_T,
      hidden_idx = hidden_positive_idx,
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
      decision_scores,
      TRUE_T
    ),

    Hidden_AUC = roc_auc(
      unlabeled_decision_scores,
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
    decision_score = decision_scores,
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
## 9. EXACT SAVED SAR-GEOMETRY VALIDATION AND FEATURE BUNDLES
## ============================================================================

same_numeric <- function(x, y, tolerance = 1e-12) {
  isTRUE(all.equal(
    as.numeric(x),
    as.numeric(y),
    tolerance = tolerance,
    check.attributes = FALSE
  ))
}

saved_sar_geometry_is_compatible <- function(geometry, seed_row) {
  if (!is.list(geometry)) return(FALSE)

  required_fields <- c(
    "settings_signature",
    "seeds",
    "T",
    "Y",
    "hidden_positive_idx",
    "observed_positive_idx",
    "positive_signal_group",
    "signal_type",
    "covariate_signal_idx",
    "network_signal_idx",
    "X_raw",
    "X_std",
    "A_network",
    "network_summary",
    "Z_network_modularity",
    "modularity_q",
    "modularity_selected_eigenvalues",
    "sar_selection_summary",
    "sar_selection_diagnostics"
  )

  if (any(!vapply(
    required_fields,
    function(name_now) !is.null(geometry[[name_now]]),
    logical(1)
  ))) {
    return(FALSE)
  }

  if (!identical(as.integer(geometry$T), TRUE_T)) return(FALSE)

  observed_y <- as.integer(geometry$Y)
  hidden_positive_idx <- sort(as.integer(geometry$hidden_positive_idx))
  observed_positive_idx <- sort(as.integer(geometry$observed_positive_idx))

  if (
    length(observed_y) != N_NODES ||
    !all(observed_y %in% c(0L, 1L)) ||
    sum(observed_y) != N_OBSERVED_POS ||
    length(hidden_positive_idx) != N_HIDDEN ||
    length(observed_positive_idx) != N_OBSERVED_POS ||
    !identical(
      sort(c(hidden_positive_idx, observed_positive_idx)),
      TRUE_POS_IDX
    ) ||
    !identical(which(observed_y == 1L), observed_positive_idx) ||
    any(observed_y[hidden_positive_idx] != 0L)
  ) {
    return(FALSE)
  }

  if (
    !identical(
      as.character(geometry$positive_signal_group),
      as.character(POSITIVE_SIGNAL_GROUP)
    ) ||
    length(geometry$signal_type) != N_NODES ||
    any(is.na(geometry$signal_type)) ||
    !identical(
      sort(as.integer(geometry$covariate_signal_idx)),
      sort(as.integer(COVARIATE_SIGNAL_IDX))
    ) ||
    !identical(
      sort(as.integer(geometry$network_signal_idx)),
      sort(as.integer(NETWORK_SIGNAL_IDX))
    )
  ) {
    return(FALSE)
  }

  required_seed_names <- c(
    "covariate_seed",
    "network_seed",
    "modularity_seed",
    "sar_seed"
  )

  if (!all(required_seed_names %in% names(geometry$seeds))) {
    return(FALSE)
  }

  saved_seeds <- as.integer(unlist(
    geometry$seeds[1L, required_seed_names, drop = FALSE]
  ))
  expected_seeds <- as.integer(unlist(
    seed_row[1L, required_seed_names, drop = FALSE]
  ))

  if (!identical(saved_seeds, expected_seeds)) return(FALSE)

  X_raw <- as.matrix(geometry$X_raw)
  X_std <- as.matrix(geometry$X_std)
  A_network <- as.matrix(geometry$A_network)
  Z_network <- as.matrix(geometry$Z_network_modularity)

  if (
    !all(dim(X_raw) == c(N_NODES, 10L)) ||
    !all(dim(X_std) == c(N_NODES, 10L)) ||
    !all(dim(A_network) == c(N_NODES, N_NODES)) ||
    nrow(Z_network) != N_NODES ||
    ncol(Z_network) < 1L ||
    as.integer(geometry$modularity_q) != ncol(Z_network) ||
    length(geometry$modularity_selected_eigenvalues) != ncol(Z_network) ||
    any(!is.finite(X_raw)) ||
    any(!is.finite(X_std)) ||
    any(!is.finite(A_network)) ||
    any(!is.finite(Z_network)) ||
    any(!is.finite(as.numeric(
      geometry$modularity_selected_eigenvalues
    )))
  ) {
    return(FALSE)
  }

  signature <- geometry$settings_signature
  if (
    !is.list(signature) ||
    is.null(signature$latent_truth) ||
    is.null(signature$label_mechanism) ||
    is.null(signature$covariate_dgp) ||
    is.null(signature$network_dgp) ||
    is.null(signature$modularity)
  ) {
    return(FALSE)
  }

  latent <- signature$latent_truth
  label_mechanism <- signature$label_mechanism
  cov_dgp <- signature$covariate_dgp
  net_dgp <- signature$network_dgp
  mod_set <- signature$modularity

  latent_ok <-
    identical(as.integer(latent$true_zero_idx), as.integer(TRUE_ZERO_IDX)) &&
    identical(as.integer(latent$true_positive_idx), as.integer(TRUE_POS_IDX)) &&
    identical(
      sort(as.integer(latent$covariates_only_idx)),
      sort(as.integer(COV_ONLY_IDX))
    ) &&
    identical(
      sort(as.integer(latent$network_only_idx)),
      sort(as.integer(NET_ONLY_IDX))
    ) &&
    identical(
      sort(as.integer(latent$both_idx)),
      sort(as.integer(BOTH_IDX))
    ) &&
    identical(
      sort(as.integer(latent$covariate_signal_idx)),
      sort(as.integer(COVARIATE_SIGNAL_IDX))
    ) &&
    identical(
      sort(as.integer(latent$network_signal_idx)),
      sort(as.integer(NETWORK_SIGNAL_IDX))
    ) &&
    as.integer(latent$n_hidden_per_replication) == N_HIDDEN &&
    as.integer(latent$n_observed_positive_per_replication) ==
      N_OBSERVED_POS

  label_ok <-
    grepl("SAR", as.character(label_mechanism$type), fixed = TRUE) &&
    as.integer(label_mechanism$exactly_n_hidden_without_replacement) ==
      N_HIDDEN &&
    same_numeric(label_mechanism$lambda, LAMBDA_SAR) &&
    same_numeric(label_mechanism$z_clip, SAR_Z_CLIP) &&
    identical(
      as.logical(label_mechanism$oracle_score_available_to_fitted_methods),
      FALSE
    )

  cov_ok <-
    same_numeric(cov_dgp$mu_linear, MU_LINEAR) &&
    same_numeric(cov_dgp$sd_linear, SD_LINEAR) &&
    same_numeric(cov_dgp$block2_x3_offset, BLOCK2_X3_OFFSET) &&
    same_numeric(cov_dgp$block2_x4_center, BLOCK2_X4_CENTER) &&
    same_numeric(cov_dgp$block2_x4_offset, BLOCK2_X4_OFFSET) &&
    same_numeric(cov_dgp$block2_zero_sd_x3, BLOCK2_ZERO_SD_X3) &&
    same_numeric(cov_dgp$block2_zero_sd_x4, BLOCK2_ZERO_SD_X4) &&
    same_numeric(cov_dgp$block2_pos_sd_x3, BLOCK2_POS_SD_X3) &&
    same_numeric(cov_dgp$block2_pos_sd_x4, BLOCK2_POS_SD_X4) &&
    same_numeric(cov_dgp$disc_radius, DISC_RADIUS) &&
    same_numeric(cov_dgp$ring_radius, RING_RADIUS) &&
    same_numeric(cov_dgp$ring_sd, RING_SD) &&
    same_numeric(cov_dgp$parabola_coef, PARABOLA_COEF) &&
    same_numeric(cov_dgp$delta_parabola, DELTA_PARABOLA) &&
    same_numeric(cov_dgp$sd_parabola_zero, SD_PARABOLA_ZERO) &&
    same_numeric(cov_dgp$sd_parabola_pos, SD_PARABOLA_POS) &&
    same_numeric(cov_dgp$rho_zero, RHO_ZERO) &&
    same_numeric(cov_dgp$rho_pos, RHO_POS)

  net_ok <-
    same_numeric(net_dgp$alpha_background, ALPHA_BACKGROUND) &&
    same_numeric(net_dgp$alpha_signal, ALPHA_FRAUD) &&
    same_numeric(net_dgp$alpha_legal, ALPHA_LEGAL) &&
    same_numeric(net_dgp$d_background, D_BACKGROUND) &&
    same_numeric(net_dgp$d_signal, D_FRAUD) &&
    same_numeric(net_dgp$d_legal, D_LEGAL) &&
    same_numeric(net_dgp$triad_closure_prob, TRIAD_CLOSURE_PROB)

  mod_ok <-
    same_numeric(mod_set$eig_tol, MOD_EIG_TOL) &&
    same_numeric(mod_set$k_max, MOD_K_MAX) &&
    (
      (is.infinite(mod_set$max_q) && is.infinite(MOD_MAX_Q)) ||
        same_numeric(mod_set$max_q, MOD_MAX_Q)
    )

  selection_summary <- as.data.frame(geometry$sar_selection_summary)
  selection_diagnostics <- as.data.frame(
    geometry$sar_selection_diagnostics
  )

  selection_ok <-
    nrow(selection_summary) == 1L &&
    nrow(selection_diagnostics) == N_TRUE_POS &&
    same_numeric(selection_summary$lambda_SAR, LAMBDA_SAR) &&
    as.integer(selection_summary$n_hidden) == N_HIDDEN &&
    as.integer(selection_summary$n_observed_positive) ==
      N_OBSERVED_POS

  isTRUE(
    latent_ok &&
      label_ok &&
      cov_ok &&
      net_ok &&
      mod_ok &&
      selection_ok
  )
}

finalize_sar_feature_bundle <- function(
    geometry,
    seed_row
) {
  X_cov <- as.matrix(geometry$X_std)
  A_network <- sanitize_network_adjacency(geometry$A_network)
  Z_net <- as.matrix(geometry$Z_network_modularity)

  observed_y <- as.integer(geometry$Y)
  hidden_positive_idx <- sort(as.integer(
    geometry$hidden_positive_idx
  ))
  observed_positive_idx <- sort(as.integer(
    geometry$observed_positive_idx
  ))
  signal_type <- as.character(geometry$signal_type)

  list(
    source = "saved_competing_SAR_geometry",
    seeds = seed_row,
    X_cov = X_cov,
    A_network = A_network,
    Z_net = Z_net,
    X_cov_net = cbind(X_cov, Z_net),
    observed_y = observed_y,
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = observed_positive_idx,
    signal_type = signal_type,
    network_summary = geometry$network_summary,
    modularity_q = as.integer(geometry$modularity_q),
    modularity_selected_eigenvalues = as.numeric(
      geometry$modularity_selected_eigenvalues
    ),
    X_raw = geometry$X_raw,
    X_center = geometry$X_center,
    X_scale = geometry$X_scale,
    covariate_profile = geometry$covariate_profile,
    sar_selection_summary = as.data.frame(
      geometry$sar_selection_summary
    ),
    sar_selection_diagnostics = as.data.frame(
      geometry$sar_selection_diagnostics
    )
  )
}

load_exact_sar_feature_bundle <- function(replicate_number) {
  seed_row <- seed_table[
    seed_table$replicate == replicate_number,
    ,
    drop = FALSE
  ]

  if (nrow(seed_row) != 1L) {
    stop("Seed lookup failed for replicate ", replicate_number, ".")
  }

  geometry_file <- file.path(
    COMPETING_GEOMETRY_DIR,
    sprintf("replicate_%03d_geometry.rds", replicate_number)
  )

  if (!file.exists(geometry_file)) {
    stop(
      "The exact competing-method SAR geometry is missing for replicate ",
      replicate_number,
      ": ",
      geometry_file,
      ". Run the SAR competing-method simulation first."
    )
  }

  geometry <- tryCatch(
    readRDS(geometry_file),
    error = function(e) NULL
  )

  if (
    is.null(geometry) ||
    !saved_sar_geometry_is_compatible(geometry, seed_row)
  ) {
    stop(
      "The saved competing-method SAR geometry failed strict validation ",
      "for replicate ",
      replicate_number,
      "."
    )
  }

  finalize_sar_feature_bundle(
    geometry = geometry,
    seed_row = seed_row
  )
}

## ============================================================================
## 10. CHECKPOINT SIGNATURE
## ============================================================================

SETTINGS_SIGNATURE <- list(
  script_version = SCRIPT_VERSION,
  master_seed = MASTER_SEED,
  n_rep = N_REP,
  n_nodes = N_NODES,

  latent_truth = list(
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
    z_clip = SAR_Z_CLIP,
    oracle_score_available_to_nnpu = FALSE
  ),

  nnpu_training_sample = list(
    class_prior = NNPU_CLASS_PRIOR,
    positive_rows = N_OBSERVED_POS,
    marginal_u_rows = N_NNPU_MARGINAL_U,
    total_training_rows = N_NNPU_TRAIN_ROWS,
    observed_positives_appear_in_P_and_U = TRUE
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
    triad_closure_prob = TRIAD_CLOSURE_PROB
  ),

  modularity = list(
    eig_tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q,
    k_max = MOD_K_MAX
  ),

  nnpu = list(
    model = "d-100-1 ReLU",
    hidden_units = NNPU_HIDDEN_UNITS,
    epochs = NNPU_EPOCHS,
    full_batch_updates = NNPU_EPOCHS,
    batch_mode = paste0(
      "30 replication-specific observed P rows and all 400 nodes ",
      "as marginal U; 430 training rows"
    ),
    surrogate_loss = "sigmoid",
    beta = NNPU_BETA,
    gamma = NNPU_GAMMA,
    optimizer = "Chainer-style Adam",
    stepsize = NNPU_STEPSIZE,
    adam_beta1 = NNPU_ADAM_BETA1,
    adam_beta2 = NNPU_ADAM_BETA2,
    adam_eps = NNPU_ADAM_EPS,
    weight_decay = NNPU_WEIGHT_DECAY,
    initialization = "LeCun normal; zero biases",
    probability_map = NNPU_PROBABILITY_MAP,
    ranking_input = paste0(
      "raw decision scores for AUC and recall ranking; ",
      "plogis(score) only for LogLoss and Brier"
    ),
    registry = method_registry
  ),

  evaluation = list(
    hidden_top_k = HIDDEN_TOP_K,
    overall_recall_k = OVERALL_RECALL_K,
    auc_and_ranking_scores = "raw nnPU decision scores",
    probability_metrics = "plogis(raw nnPU decision score)",
    predictive_interval_levels = PREDICTIVE_INTERVAL_LEVELS,
    uncertainty = "NA for both nnPU versions"
  )
)

## ============================================================================
## 11. ONE COMPLETE REPLICATION
## ============================================================================

run_one_replication <- function(replicate_number) {
  replicate_number <- as.integer(replicate_number)

  seed_row <- seed_table[
    seed_table$replicate == replicate_number,
    ,
    drop = FALSE
  ]

  if (nrow(seed_row) != 1L) {
    stop("Seed-table lookup failed for replicate ", replicate_number, ".")
  }

  checkpoint_file <- file.path(
    CHECKPOINT_DIR,
    sprintf("replicate_%03d.rds", replicate_number)
  )

  checkpoint <- NULL

  if (
    isTRUE(RESUME_REPLICATES) &&
    !isTRUE(FORCE_FRESH_RUN) &&
    file.exists(checkpoint_file)
  ) {
    candidate <- tryCatch(
      readRDS(checkpoint_file),
      error = function(e) NULL
    )

    if (
      is.list(candidate) &&
      identical(candidate$settings_signature, SETTINGS_SIGNATURE) &&
      identical(as.integer(candidate$replicate), replicate_number)
    ) {
      checkpoint <- candidate
    }
  }

  if (!is.null(checkpoint) && isTRUE(checkpoint$complete)) {
    return(checkpoint)
  }

  replication_start <- Sys.time()

  bundle <- load_exact_sar_feature_bundle(
    replicate_number = replicate_number
  )

  X_cov <- bundle$X_cov
  Z_net <- bundle$Z_net
  X_cov_net <- bundle$X_cov_net
  observed_y <- bundle$observed_y
  hidden_positive_idx <- bundle$hidden_positive_idx
  observed_positive_idx <- bundle$observed_positive_idx
  signal_type <- bundle$signal_type

  selection_summary_now <- bundle$sar_selection_summary
  selection_diagnostics_now <- bundle$sar_selection_diagnostics

  if (!"replicate" %in% names(selection_summary_now)) {
    selection_summary_now$replicate <- replicate_number
  }
  if (!"replicate" %in% names(selection_diagnostics_now)) {
    selection_diagnostics_now$replicate <- replicate_number
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
      sar_selection_summary = NULL,
      sar_selection_diagnostics = NULL,
      label_design = NULL,
      time_start = replication_start,
      time_end = NULL,
      runtime_sec = NA_real_
    )
  }

  checkpoint$geometry_summary <- cbind(
    data.frame(
      scenario = SCENARIO_NAME,
      replicate = replicate_number,
      geometry_source = bundle$source,
      covariate_seed = seed_row$covariate_seed,
      network_seed = seed_row$network_seed,
      modularity_seed = seed_row$modularity_seed,
      sar_seed = seed_row$sar_seed,
      lambda_SAR = LAMBDA_SAR,
      n = N_NODES,
      n_true_zero = N_TRUE_ZERO,
      n_true_positive = N_TRUE_POS,
      n_observed_positive = length(observed_positive_idx),
      n_hidden_positive = length(hidden_positive_idx),
      hidden_covariates_only = sum(
        hidden_positive_idx %in% COV_ONLY_IDX
      ),
      hidden_network_only = sum(
        hidden_positive_idx %in% NET_ONLY_IDX
      ),
      hidden_both = sum(
        hidden_positive_idx %in% BOTH_IDX
      ),
      n_covariate_signal = length(COVARIATE_SIGNAL_IDX),
      n_network_signal = length(NETWORK_SIGNAL_IDX),
      n_both_signal = length(intersect(
        COVARIATE_SIGNAL_IDX,
        NETWORK_SIGNAL_IDX
      )),
      stringsAsFactors = FALSE
    ),
    bundle$network_summary
  )

  checkpoint$modularity_summary <- data.frame(
    scenario = SCENARIO_NAME,
    replicate = replicate_number,
    geometry_source = bundle$source,
    q_network = bundle$modularity_q,
    selected_eigenvalues = paste(
      signif(bundle$modularity_selected_eigenvalues, 8),
      collapse = ";"
    ),
    n_covariates = ncol(X_cov),
    n_network_features = ncol(Z_net),
    n_covariates_plus_network = ncol(X_cov_net),
    stringsAsFactors = FALSE
  )

  checkpoint$sar_selection_summary <- selection_summary_now
  checkpoint$sar_selection_diagnostics <- selection_diagnostics_now
  checkpoint$label_design <- list(
    observed_y = observed_y,
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = observed_positive_idx,
    signal_type = signal_type
  )

  if (isTRUE(SAVE_GENERATED_GEOMETRIES)) {
    saveRDS(
      list(
        settings_signature = SETTINGS_SIGNATURE,
        source = bundle$source,
        scenario = SCENARIO_NAME,
        replicate = replicate_number,
        seeds = seed_row,
        T = TRUE_T,
        Y = observed_y,
        hidden_positive_idx = hidden_positive_idx,
        observed_positive_idx = observed_positive_idx,
        positive_signal_group = POSITIVE_SIGNAL_GROUP,
        signal_type = signal_type,
        covariate_signal_idx = COVARIATE_SIGNAL_IDX,
        network_signal_idx = NETWORK_SIGNAL_IDX,
        covariate_profile = bundle$covariate_profile,
        X_raw = bundle$X_raw,
        X_std = X_cov,
        X_center = bundle$X_center,
        X_scale = bundle$X_scale,
        A_network = bundle$A_network,
        network_summary = bundle$network_summary,
        Z_network_modularity = Z_net,
        modularity_q = bundle$modularity_q,
        modularity_selected_eigenvalues =
          bundle$modularity_selected_eigenvalues,
        sar_selection_summary = selection_summary_now,
        sar_selection_diagnostics = selection_diagnostics_now
      ),
      file.path(
        GEOMETRY_DIR,
        sprintf("replicate_%03d_geometry.rds", replicate_number)
      )
    )
  }

  cat(sprintf(
    paste0(
      "[%s] Replicate %03d/%03d | SAR hidden composition: ",
      "cov=%d, net=%d, both=%d\n"
    ),
    timestamp_now(),
    replicate_number,
    N_REP,
    sum(hidden_positive_idx %in% COV_ONLY_IDX),
    sum(hidden_positive_idx %in% NET_ONLY_IDX),
    sum(hidden_positive_idx %in% BOTH_IDX)
  ))

  for (method_index in seq_len(nrow(method_registry))) {
    method_id <- method_registry$method_id[method_index]
    method_label <- method_registry$method[method_index]

    if (method_id %in% names(checkpoint$method_results)) {
      cat(sprintf(
        "[%s] Replicate %03d/%03d | reusing %s\n",
        timestamp_now(),
        replicate_number,
        N_REP,
        method_label
      ))
      next
    }

    cat(sprintf(
      "[%s] Replicate %03d/%03d | fitting %s\n",
      timestamp_now(),
      replicate_number,
      N_REP,
      method_label
    ))

    method_start <- Sys.time()
    method_error <- NA_character_

    method_output <- tryCatch(
      run_nnpu_method(
        method_id = method_id,
        X_cov = X_cov,
        X_cov_net = X_cov_net,
        y = observed_y,
        replicate_seed_row = seed_row
      ),
      error = function(e) {
        method_error <<- conditionMessage(e)
        NULL
      }
    )

    method_runtime_sec <- as.numeric(difftime(
      Sys.time(),
      method_start,
      units = "secs"
    ))

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

    evaluation_now <- evaluate_method_scores(
      method_output = method_output,
      method_id = method_id,
      method_label = method_label,
      replicate_number = replicate_number,
      observed_y = observed_y,
      hidden_positive_idx = hidden_positive_idx,
      signal_type = signal_type
    )

    evaluation_now$metrics$runtime_sec <- method_runtime_sec

    evaluation_now$training_trace <- method_output$training_trace
    evaluation_now$training_trace$scenario <- SCENARIO_NAME
    evaluation_now$training_trace$replicate <- replicate_number
    evaluation_now$training_trace$method_id <- method_id
    evaluation_now$training_trace$method <- method_label
    evaluation_now$training_trace <- evaluation_now$training_trace[
      ,
      c(
        "scenario",
        "replicate",
        "method_id",
        "method",
        setdiff(
          names(evaluation_now$training_trace),
          c("scenario", "replicate", "method_id", "method")
        )
      ),
      drop = FALSE
    ]

    evaluation_now$training_diagnostics <-
      method_output$training_diagnostics
    evaluation_now$training_diagnostics$scenario <- SCENARIO_NAME
    evaluation_now$training_diagnostics$replicate <- replicate_number
    evaluation_now$training_diagnostics$method_id <- method_id
    evaluation_now$training_diagnostics$method <- method_label
    evaluation_now$training_diagnostics$runtime_sec <-
      method_runtime_sec
    evaluation_now$training_diagnostics$hidden_covariates_only <-
      sum(hidden_positive_idx %in% COV_ONLY_IDX)
    evaluation_now$training_diagnostics$hidden_network_only <-
      sum(hidden_positive_idx %in% NET_ONLY_IDX)
    evaluation_now$training_diagnostics$hidden_both <-
      sum(hidden_positive_idx %in% BOTH_IDX)

    evaluation_now$training_diagnostics <-
      evaluation_now$training_diagnostics[
        ,
        c(
          "scenario",
          "replicate",
          "method_id",
          "method",
          setdiff(
            names(evaluation_now$training_diagnostics),
            c("scenario", "replicate", "method_id", "method")
          )
        ),
        drop = FALSE
      ]

    checkpoint$method_results[[method_id]] <- evaluation_now
    checkpoint$method_failures[[method_id]] <- NULL
    saveRDS(checkpoint, checkpoint_file)

    rm(method_output, evaluation_now)
    invisible(gc())
  }

  missing_method_ids <- setdiff(
    method_registry$method_id,
    names(checkpoint$method_results)
  )

  if (length(missing_method_ids) > 0L) {
    checkpoint$complete <- FALSE
    saveRDS(checkpoint, checkpoint_file)
    stop(
      "Replication ",
      replicate_number,
      " did not complete every nnPU method. Missing: ",
      paste(missing_method_ids, collapse = ", "),
      "."
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

  rm(bundle, X_cov, Z_net, X_cov_net)
  invisible(gc())

  checkpoint
}

## ============================================================================
## 12. RUN ALL 100 REPLICATIONS
## ============================================================================

replication_results <- vector("list", N_REP)
simulation_start <- Sys.time()
n_loaded_complete <- 0L
n_computed_now <- 0L

for (replicate_number in seq_len(N_REP)) {
  checkpoint_file <- file.path(
    CHECKPOINT_DIR,
    sprintf("replicate_%03d.rds", replicate_number)
  )
  was_complete_before <- FALSE

  if (isTRUE(RESUME_REPLICATES) && !isTRUE(FORCE_FRESH_RUN) &&
      file.exists(checkpoint_file)) {
    prior_checkpoint <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
    was_complete_before <-
      is.list(prior_checkpoint) &&
      identical(prior_checkpoint$settings_signature, SETTINGS_SIGNATURE) &&
      isTRUE(prior_checkpoint$complete)
  }

  replication_results[[replicate_number]] <-
    run_one_replication(replicate_number)

  if (was_complete_before) {
    n_loaded_complete <- n_loaded_complete + 1L
  } else {
    n_computed_now <- n_computed_now + 1L
  }

  elapsed_sec <- as.numeric(difftime(
    Sys.time(), simulation_start, units = "secs"
  ))
  projected_total_sec <- elapsed_sec / replicate_number * N_REP

  cat(sprintf(
    "[%s] Completed replication %03d/%03d | elapsed %s | projected total %s\n",
    timestamp_now(),
    replicate_number,
    N_REP,
    format_duration(elapsed_sec),
    format_duration(projected_total_sec)
  ))
}

simulation_end <- Sys.time()

## ============================================================================
## 13. COMBINE RESULTS
## ============================================================================

if (length(replication_results) != N_REP ||
    any(!vapply(
      replication_results,
      function(x) is.list(x) && isTRUE(x$complete),
      logical(1)
    ))) {
  stop("Not all 100 nnPU simulation replications completed successfully.")
}

metrics_by_replicate <- do.call(
  rbind,
  lapply(replication_results, function(replication_result) {
    do.call(rbind, lapply(method_registry$method_id, function(method_id) {
      replication_result$method_results[[method_id]]$metrics
    }))
  })
)
row.names(metrics_by_replicate) <- NULL

scores_all_nodes <- do.call(
  rbind,
  lapply(replication_results, function(replication_result) {
    do.call(rbind, lapply(method_registry$method_id, function(method_id) {
      replication_result$method_results[[method_id]]$node_scores
    }))
  })
)
row.names(scores_all_nodes) <- NULL

training_diagnostics <- do.call(
  rbind,
  lapply(replication_results, function(replication_result) {
    do.call(rbind, lapply(method_registry$method_id, function(method_id) {
      replication_result$method_results[[method_id]]$training_diagnostics
    }))
  })
)
row.names(training_diagnostics) <- NULL

training_trace <- if (isTRUE(SAVE_TRAINING_TRACES)) {
  out <- do.call(
    rbind,
    lapply(replication_results, function(replication_result) {
      do.call(rbind, lapply(method_registry$method_id, function(method_id) {
        replication_result$method_results[[method_id]]$training_trace
      }))
    })
  )
  row.names(out) <- NULL
  out
} else {
  NULL
}

geometry_summary <- do.call(
  rbind,
  lapply(replication_results, `[[`, "geometry_summary")
)
row.names(geometry_summary) <- NULL

modularity_summary <- do.call(
  rbind,
  lapply(replication_results, `[[`, "modularity_summary")
)
row.names(modularity_summary) <- NULL


sar_selection_summary <- do.call(
  rbind,
  lapply(replication_results, `[[`, "sar_selection_summary")
)
row.names(sar_selection_summary) <- NULL

sar_selection_diagnostics <- do.call(
  rbind,
  lapply(replication_results, `[[`, "sar_selection_diagnostics")
)
row.names(sar_selection_diagnostics) <- NULL

if (
  nrow(sar_selection_summary) != N_REP ||
  nrow(sar_selection_diagnostics) != N_REP * N_TRUE_POS
) {
  stop("The combined SAR-selection diagnostics have invalid dimensions.")
}

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

failure_log_parts <- Filter(
  Negate(is.null),
  lapply(replication_results, function(replication_result) {
    failures_now <- replication_result$method_failures
    if (length(failures_now) == 0L) return(NULL)
    do.call(rbind, failures_now)
  })
)

failure_log <- if (length(failure_log_parts) == 0L) {
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
  do.call(rbind, failure_log_parts)
}

if (nrow(metrics_by_replicate) != N_REP * nrow(method_registry)) {
  stop("The metric table does not contain 100 x 2 rows.")
}
if (nrow(scores_all_nodes) != N_REP * nrow(method_registry) * N_NODES) {
  stop("The node-score table does not contain 100 x 2 x 400 rows.")
}

method_order <- method_registry$method
metrics_by_replicate$method <- factor(
  metrics_by_replicate$method,
  levels = method_order
)
metrics_by_replicate <- metrics_by_replicate[
  order(metrics_by_replicate$method, metrics_by_replicate$replicate),
  ,
  drop = FALSE
]
metrics_by_replicate$method <- as.character(metrics_by_replicate$method)

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
scores_all_nodes$method <- as.character(scores_all_nodes$method)

write.csv(
  metrics_by_replicate,
  file.path(OUT_DIR, "nnpu_metrics_by_replicate.csv"),
  row.names = FALSE
)
write.csv(
  scores_all_nodes,
  file.path(OUT_DIR, "nnpu_all_node_probabilities.csv"),
  row.names = FALSE
)
write.csv(
  training_diagnostics,
  file.path(OUT_DIR, "nnpu_training_diagnostics.csv"),
  row.names = FALSE
)
if (!is.null(training_trace)) {
  write.csv(
    training_trace,
    file.path(OUT_DIR, "nnpu_training_trace_all_epochs.csv"),
    row.names = FALSE
  )
}
write.csv(
  geometry_summary,
  file.path(OUT_DIR, "nnpu_geometry_summary.csv"),
  row.names = FALSE
)
write.csv(
  modularity_summary,
  file.path(OUT_DIR, "nnpu_modularity_summary.csv"),
  row.names = FALSE
)
write.csv(
  sar_selection_summary,
  file.path(OUT_DIR, "nnpu_SAR_selection_summary_by_replicate.csv"),
  row.names = FALSE
)
write.csv(
  sar_selection_diagnostics,
  file.path(OUT_DIR, "nnpu_SAR_selection_diagnostics_all_positives.csv"),
  row.names = FALSE
)
write.csv(
  sar_selection_mc_summary,
  file.path(OUT_DIR, "nnpu_SAR_selection_Monte_Carlo_summary.csv"),
  row.names = FALSE
)
write.csv(
  failure_log,
  file.path(OUT_DIR, "nnpu_failure_log.csv"),
  row.names = FALSE
)

## ============================================================================
## 14. MONTE CARLO SUMMARIES AND PAPER TABLES
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

format_mean_sd <- function(mean_value, sd_value, digits = 3L) {
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

build_metric_summaries <- function(
    metric_data,
    method_order_now,
    output_directory,
    file_prefix
) {
  summary_rows <- list()
  summary_index <- 1L

  for (method_name in method_order_now) {
    method_data <- metric_data[
      metric_data$method == method_name,
      ,
      drop = FALSE
    ]

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

  main_table_means <- data.frame(
    scenario = SCENARIO_NAME,
    method = method_order_now,
    stringsAsFactors = FALSE
  )

  paper_table_mean_sd <- data.frame(
    Method = method_order_now,
    stringsAsFactors = FALSE
  )

  for (metric_name in main_metric_names) {
    metric_rows <- metric_summary_long[
      metric_summary_long$metric == metric_name,
      ,
      drop = FALSE
    ]
    metric_rows <- metric_rows[
      match(method_order_now, metric_rows$method),
      ,
      drop = FALSE
    ]

    main_table_means[[metric_name]] <- metric_rows$mean
    paper_table_mean_sd[[paper_column_labels[metric_name]]] <- mapply(
      format_mean_sd,
      metric_rows$mean,
      metric_rows$sd,
      MoreArgs = list(digits = 3L),
      USE.NAMES = FALSE
    )
  }

  write.csv(
    metric_summary_long,
    file.path(
      output_directory,
      paste0(file_prefix, "_metric_summary_long.csv")
    ),
    row.names = FALSE
  )
  write.csv(
    main_table_means,
    file.path(
      output_directory,
      paste0(file_prefix, "_FINAL_MAIN_TABLE_means.csv")
    ),
    row.names = FALSE
  )
  write.csv(
    paper_table_mean_sd,
    file.path(
      output_directory,
      paste0(file_prefix, "_FINAL_PAPER_TABLE_mean_SD.csv")
    ),
    row.names = FALSE,
    na = ""
  )

  list(
    metric_summary_long = metric_summary_long,
    main_table_means = main_table_means,
    paper_table_mean_sd = paper_table_mean_sd
  )
}

nnpu_summary_outputs <- build_metric_summaries(
  metric_data = metrics_by_replicate,
  method_order_now = method_order,
  output_directory = OUT_DIR,
  file_prefix = "nnpu"
)

metric_summary_long <- nnpu_summary_outputs$metric_summary_long
main_table_means <- nnpu_summary_outputs$main_table_means
paper_table_mean_sd <- nnpu_summary_outputs$paper_table_mean_sd

## Optional combined seven-original-plus-two-nnPU table.
combined_summary_outputs <- NULL
competitor_metrics_file <- file.path(
  COMPETING_OUT_DIR,
  "competing_methods_metrics_by_replicate.csv"
)

if (isTRUE(MERGE_WITH_EXISTING_COMPETITOR_RESULTS) &&
    file.exists(competitor_metrics_file)) {
  competitor_metrics <- tryCatch(
    read.csv(competitor_metrics_file, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  required_metric_columns <- names(metrics_by_replicate)
  if (!is.null(competitor_metrics) &&
      all(required_metric_columns %in% names(competitor_metrics))) {
    competitor_metrics <- competitor_metrics[
      ,
      required_metric_columns,
      drop = FALSE
    ]
    combined_metrics <- rbind(
      competitor_metrics,
      metrics_by_replicate
    )
    combined_method_order <- c(
      unique(competitor_metrics$method),
      method_order
    )
    combined_method_order <- unique(combined_method_order)

    write.csv(
      combined_metrics,
      file.path(
        OUT_DIR,
        "all_competitors_plus_nnpu_metrics_by_replicate.csv"
      ),
      row.names = FALSE
    )

    combined_summary_outputs <- build_metric_summaries(
      metric_data = combined_metrics,
      method_order_now = combined_method_order,
      output_directory = OUT_DIR,
      file_prefix = "all_competitors_plus_nnpu"
    )
  } else {
    warning(
      "The existing competing-method metric file was found but has an ",
      "incompatible schema; the combined table was skipped."
    )
  }
}

## ============================================================================
## 15. SAVE COMPLETE OUTPUT BUNDLE
## ============================================================================

total_runtime_sec <- as.numeric(difftime(
  simulation_end,
  simulation_start,
  units = "secs"
))

run_config <- list(
  script_version = SCRIPT_VERSION,
  created_at = timestamp_now(),
  scenario = SCENARIO_NAME,
  settings_signature = SETTINGS_SIGNATURE,
  output_directory = OUT_DIR,
  competing_geometry_directory = COMPETING_GEOMETRY_DIR,
  require_saved_competing_geometries =
    REQUIRE_SAVED_COMPETING_GEOMETRIES,
  n_rep = N_REP,
  n_nodes = N_NODES,
  n_evaluation_unlabeled = N_UNLABELED,
  n_nnpu_positive_rows = N_OBSERVED_POS,
  n_nnpu_marginal_u_rows = N_NNPU_MARGINAL_U,
  n_nnpu_training_rows = N_NNPU_TRAIN_ROWS,
  n_full_batch_updates = NNPU_EPOCHS,
  class_prior = NNPU_CLASS_PRIOR,

  latent_truth = list(
    true_labels = TRUE_T,
    true_zero_idx = TRUE_ZERO_IDX,
    true_positive_idx = TRUE_POS_IDX,
    positive_signal_group = POSITIVE_SIGNAL_GROUP,
    covariates_only_idx = COV_ONLY_IDX,
    network_only_idx = NET_ONLY_IDX,
    both_idx = BOTH_IDX,
    covariate_signal_idx = COVARIATE_SIGNAL_IDX,
    network_signal_idx = NETWORK_SIGNAL_IDX
  ),

  label_mechanism = list(
    type = "covariate-driven SAR",
    lambda = LAMBDA_SAR,
    z_clip = SAR_Z_CLIP,
    n_hidden_per_replication = N_HIDDEN,
    n_observed_positive_per_replication = N_OBSERVED_POS,
    oracle_score_used_for_nnpu = FALSE
  ),

  sar_selection_summary = sar_selection_summary,
  sar_selection_Monte_Carlo_summary = sar_selection_mc_summary,

  methods = method_registry,

  nnpu_probability_note = paste0(
    "The official nnPU classifier returns decision scores. AUC and recall ",
    "ranking use raw scores. LogLoss and Brier use plogis(score), without ",
    "post-hoc calibration or access to TRUE_T."
  ),

  nnpu_training_sample_note = paste0(
    "In every SAR replication, each update uses the 30 replication-specific ",
    "observed positives as P and all 400 nodes as marginal U, for 430 rows. ",
    "The observed positives appear once as P and once inside U."
  ),

  predictive_uncertainty_note = paste0(
    "No model-based predictive distribution is available in the official ",
    "point-estimation nnPU implementation; all coverage and interval-score ",
    "entries are NA."
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
  file.path(OUT_DIR, "nnpu_SAR_simulation_run_config.rds")
)

save(
  replication_results,
  metrics_by_replicate,
  scores_all_nodes,
  training_diagnostics,
  training_trace,
  geometry_summary,
  modularity_summary,
  sar_selection_summary,
  sar_selection_diagnostics,
  sar_selection_mc_summary,
  failure_log,
  metric_summary_long,
  main_table_means,
  paper_table_mean_sd,
  combined_summary_outputs,
  method_registry,
  seed_table,
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
  run_config,
  SETTINGS_SIGNATURE,
  file = file.path(
    OUT_DIR,
    "SIMULATION_100_NNPU_SAR_COMPLETE_RESULTS.RData"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo_simulation_nnpu_SAR.txt")
)

## ============================================================================
## 16. FINAL CONSOLE OUTPUT
## ============================================================================

cat("\n================ RUN ACCOUNTING ================\n")
cat("Requested replications :", N_REP, "\n")
cat("Computed in this run   :", n_computed_now, "\n")
cat("Loaded from checkpoint :", n_loaded_complete, "\n")
cat("Total runtime          :", format_duration(total_runtime_sec), "\n")

cat("\n================ SAR SELECTION MONTE CARLO SUMMARY ================\n")
print(sar_selection_mc_summary, row.names = FALSE)

cat("\n================ FINAL nnPU SAR MONTE CARLO MEANS ================\n")
print(main_table_means, row.names = FALSE)

cat("\n================ FINAL nnPU SAR PAPER TABLE: MEAN (SD) ================\n")
print(paper_table_mean_sd, row.names = FALSE)

if (!is.null(combined_summary_outputs)) {
  cat("\n================ SAR COMPETITORS + nnPU: MEAN (SD) ================\n")
  print(combined_summary_outputs$paper_table_mean_sd, row.names = FALSE)
}

cat("\nAll outputs were written to:\n")
cat(normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE), "\n")
cat("\n================ DONE ================\n")
