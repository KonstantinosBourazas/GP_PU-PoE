
## ============================================================================
## SIMULATION STUDY: ALL COMPETING METHODS ON 100 MIXED-SIGNAL SAR DATA SETS
## ----------------------------------------------------------------------------
## This script is the SAR counterpart of the SCAR competing-method script.
## Every component of the data-generating process, every fitted method, every
## metric, every tuning constant and every output table is retained unchanged.
##
## The only substantive change is the observed-label mechanism:
##
##   - the same 400-node mixed-signal truth is generated in every replication;
##   - there are still 45 true positives, 30 observed positives and 15 hidden
##     positives;
##   - the 15 hidden positives are selected without replacement under a
##     covariate-driven SAR mechanism with lambda = 2;
##   - positives with smaller oracle covariate scores, i.e. those whose
##     covariates look more like the legal/background mechanism, receive larger
##     hiding weights;
##   - the oracle score is used only to generate missing labels and is never
##     supplied to any fitted method.
##
## Mixed-signal truth retained from the SCAR script
## ------------------------------------------------
##   True zeros: 1:355
##
##   Among the 45 true positives:
##     15 carry covariate signal only;
##     15 carry network signal only;
##     15 carry both covariate and network signal.
##
## The exact signal-node sets are unchanged from the SCAR script. Only the
## hidden/observed status is regenerated in each replication through SAR.
##
## Competing methods retained unchanged
## ------------------------------------
##   Ada, Ada + Net, WLR, WLR + Net, INLA, INLA + Net, GNN.
##
## Metrics retained unchanged
## --------------------------
##   Overall and hidden AUC;
##   overall and hidden LogLoss;
##   overall and hidden Brier score;
##   Overall Recall@45 and Hidden Recall@15;
##   80% and 95% predictive coverage and interval scores where available.
##
## Reproducibility and restart safety
## ----------------------------------
## Covariate, network, modularity, SAR-selection, fitting and uncertainty seeds
## are deterministic. A new script version and output directory prevent reuse
## of checkpoints from the SCAR study.
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


SCRIPT_VERSION <- "simulation_competitors_mixed_signal_SAR_lambda2_all_metrics_v1_0"

MASTER_SEED <- 20260718L
N_REP <- 100L

N_NODES <- 400L

## --------------------------------------------------------------------------
## 1A. Fixed latent truth and fixed mixed-signal geometry
## --------------------------------------------------------------------------

TRUE_ZERO_IDX <- 1:355
TRUE_POS_IDX <- 356:400

## These three sets are the same signal sets implied by the original
## SCAR script. Their hidden/observed status is no longer fixed by index.
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

## --------------------------------------------------------------------------
## 1B. SAR label mechanism
## --------------------------------------------------------------------------

LAMBDA_SAR <- 2.00
SAR_Z_CLIP <- 3.00

SCENARIO_NAME <- paste0(
  "Fixed mixed-signal covariate-driven SAR design, lambda = ",
  format(LAMBDA_SAR, trim = TRUE)
)

OUT_DIR <- file.path(
  getwd(),
  "SIMULATION_100_COMPETING_METHODS_MIXED_SIGNAL_SAR_LAMBDA_2"
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

## Safe restart behavior.
RESUME_REPLICATES <- TRUE
FORCE_FRESH_RUN <- FALSE
FAIL_IF_ANY_METHOD_FAILS <- TRUE

## Save X, A, labels and signal metadata for all replications.
SAVE_GENERATED_GEOMETRIES <- TRUE

## --------------------------------------------------------------------------
## 1C. Evaluation settings
## --------------------------------------------------------------------------

HIDDEN_TOP_K <- 15L
OVERALL_RECALL_K <- 45L

PREDICTIVE_INTERVAL_LEVELS <- c(0.80, 0.95)
N_PREDICTIVE_DRAWS <- 2000L

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
## 1F. Network modularity-feature settings
## --------------------------------------------------------------------------

MOD_EIG_TOL <- 1e-8
MOD_MAX_Q <- Inf
MOD_K_MAX <- 20L

## --------------------------------------------------------------------------
## 1G. Competing-method settings
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
## 1H. Reproducibility settings
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
  ada_seed = as.integer(
    MASTER_SEED + 400000L + 1000L * seq_len(N_REP) + 101L
  ),
  wlr_seed = as.integer(
    MASTER_SEED + 500000L + 1000L * seq_len(N_REP) + 103L
  ),
  inla_seed = as.integer(
    MASTER_SEED + 600000L + 1000L * seq_len(N_REP) + 107L
  ),
  gnn_seed = as.integer(
    MASTER_SEED + 700000L + 1000L * seq_len(N_REP) + 109L
  ),
  wlr_predictive_seed = as.integer(
    MASTER_SEED + 800000L + 1000L * seq_len(N_REP) + 113L
  ),
  wlr_net_predictive_seed = as.integer(
    MASTER_SEED + 900000L + 1000L * seq_len(N_REP) + 127L
  ),
  inla_predictive_seed = as.integer(
    MASTER_SEED + 1000000L + 1000L * seq_len(N_REP) + 131L
  ),
  inla_net_predictive_seed = as.integer(
    MASTER_SEED + 1100000L + 1000L * seq_len(N_REP) + 137L
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
      seed_table[, seed_columns, drop = FALSE]
    ) > .Machine$integer.max
  )
) {
  stop("At least one deterministic seed exceeds the R integer range.")
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
  stop("The requested simulation study requires N_REP = 100.")
}

if (N_NODES != 400L) {
  stop("The requested simulation study requires exactly 400 observations.")
}

if (!identical(TRUE_ZERO_IDX, 1:355)) {
  stop("TRUE_ZERO_IDX must be exactly 1:355.")
}

if (!identical(TRUE_POS_IDX, 356:400)) {
  stop("TRUE_POS_IDX must be exactly 356:400.")
}

if (
  length(COV_ONLY_IDX) != 15L ||
  length(NET_ONLY_IDX) != 15L ||
  length(BOTH_IDX) != 15L
) {
  stop("Each mixed-signal positive group must contain exactly 15 nodes.")
}

if (
  length(unique(c(
    COV_ONLY_IDX,
    NET_ONLY_IDX,
    BOTH_IDX
  ))) != N_TRUE_POS ||
  !identical(
    sort(c(
      COV_ONLY_IDX,
      NET_ONLY_IDX,
      BOTH_IDX
    )),
    TRUE_POS_IDX
  )
) {
  stop("The three positive signal groups must partition TRUE_POS_IDX.")
}

if (!identical(
  sort(COVARIATE_SIGNAL_IDX),
  sort(c(
    356:360,
    366:380,
    391:400
  ))
)) {
  stop(
    "COVARIATE_SIGNAL_IDX no longer matches the original SCAR geometry."
  )
}

if (!identical(
  sort(NETWORK_SIGNAL_IDX),
  sort(c(
    361:370,
    381:400
  ))
)) {
  stop(
    "NETWORK_SIGNAL_IDX no longer matches the original SCAR geometry."
  )
}

if (N_TRUE_POS != 45L) {
  stop("There must be exactly 45 true positives.")
}

if (N_HIDDEN != 15L) {
  stop("There must be exactly 15 hidden positives per replication.")
}

if (N_OBSERVED_POS != 30L) {
  stop("There must be exactly 30 observed positives per replication.")
}

if (N_TRUE_ZERO != 355L) {
  stop("There must be exactly 355 true zeros.")
}

if (N_UNLABELED != 370L) {
  stop("There must be exactly 370 unlabeled observations.")
}

if (length(COVARIATE_SIGNAL_IDX) != 30L) {
  stop("There must be exactly 30 covariate-signal observations.")
}

if (length(NETWORK_SIGNAL_IDX) != 30L) {
  stop("There must be exactly 30 network-signal observations.")
}

if (length(intersect(
  COVARIATE_SIGNAL_IDX,
  NETWORK_SIGNAL_IDX
)) != 15L) {
  stop("Exactly 15 observations must carry both signals.")
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
  !isTRUE(all.equal(PREDICTIVE_INTERVAL_LEVELS, c(0.80, 0.95)))
) {
  stop("PREDICTIVE_INTERVAL_LEVELS must equal c(0.80, 0.95).")
}

if (
  any(!is.finite(PREDICTIVE_INTERVAL_LEVELS)) ||
  any(PREDICTIVE_INTERVAL_LEVELS <= 0) ||
  any(PREDICTIVE_INTERVAL_LEVELS >= 1)
) {
  stop("Every predictive interval level must lie strictly between 0 and 1.")
}

if (
  length(N_PREDICTIVE_DRAWS) != 1L ||
  !is.finite(N_PREDICTIVE_DRAWS) ||
  N_PREDICTIVE_DRAWS < 1000L
) {
  stop("N_PREDICTIVE_DRAWS must be one integer of at least 1000.")
}

N_PREDICTIVE_DRAWS <- as.integer(N_PREDICTIVE_DRAWS)

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

cat("\n================ FIXED LATENT TRUTH AND SAR LABEL DESIGN ================\n")
cat("Replications              :", N_REP, "\n")
cat("Observations              :", N_NODES, "\n")
cat("True zeros               :", N_TRUE_ZERO, "\n")
cat("True positives           :", N_TRUE_POS, "\n")
cat("Observed positives       :", N_OBSERVED_POS, "per replication\n")
cat("Hidden positives         :", N_HIDDEN, "per replication\n")
cat("Unlabeled observations   :", N_UNLABELED, "per replication\n")
cat("Covariate-signal nodes   :", length(COVARIATE_SIGNAL_IDX), "\n")
cat("Network-signal nodes     :", length(NETWORK_SIGNAL_IDX), "\n")
cat(
  "Both-signal nodes        :",
  length(intersect(
    COVARIATE_SIGNAL_IDX,
    NETWORK_SIGNAL_IDX
  )),
  "\n"
)
cat("SAR lambda               :", LAMBDA_SAR, "\n")
cat("SAR z-score clipping     :", SAR_Z_CLIP, "\n")
cat(
  "Predictive interval levels:",
  paste(PREDICTIVE_INTERVAL_LEVELS, collapse = ", "),
  "\n"
)
cat("Predictive draws/method  :", N_PREDICTIVE_DRAWS, "\n")


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
## 4. COVARIATE-DRIVEN SAR LABEL MECHANISM
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
    replicate_number,
    observed_y,
    hidden_positive_idx,
    signal_type
) {
  if (!is.list(method_output) || is.null(method_output$probabilities)) {
    stop(method_label, " did not return the required method-output list.")
  }

  observed_y <- as.integer(observed_y)
  hidden_positive_idx <- sort(
    unique(
      as.integer(hidden_positive_idx)
    )
  )

  if (
    length(observed_y) != N_NODES ||
    !all(observed_y %in% c(0L, 1L)) ||
    sum(observed_y) != N_OBSERVED_POS
  ) {
    stop("The replication-specific observed-label vector is invalid.")
  }

  if (
    length(hidden_positive_idx) != N_HIDDEN ||
    any(!hidden_positive_idx %in% TRUE_POS_IDX) ||
    any(observed_y[hidden_positive_idx] != 0L)
  ) {
    stop("The replication-specific hidden-positive indices are invalid.")
  }

  if (
    length(signal_type) != N_NODES ||
    any(is.na(signal_type))
  ) {
    stop("The replication-specific signal_type vector is invalid.")
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

  unlabeled_idx <- which(
    observed_y == 0L
  )

  hidden_truth_unlabeled <- as.integer(
    TRUE_T[unlabeled_idx] == 1L
  )

  if (
    length(unlabeled_idx) != N_UNLABELED ||
    sum(hidden_truth_unlabeled) != N_HIDDEN ||
    !identical(
      sort(
        unlabeled_idx[
          hidden_truth_unlabeled == 1L
        ]
      ),
      hidden_positive_idx
    )
  ) {
    stop("The unlabeled SAR evaluation set is inconsistent.")
  }

  overall_rank <- rank(
    -probabilities,
    ties.method = "average"
  )

  unlabeled_probabilities <- probabilities[
    unlabeled_idx
  ]

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
    predictive_80 <- empty_predictive_label_metrics(
      N_NODES
    )

    predictive_95 <- empty_predictive_label_metrics(
      N_NODES
    )
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
      probabilities,
      TRUE_T
    ),

    Hidden_AUC = roc_auc(
      unlabeled_probabilities,
      hidden_truth_unlabeled
    ),

    Hidden_at_15 = hidden_at_15,

    Hidden_Recall_at_15 =
      hidden_at_15 /
      N_HIDDEN,

    Overall_Recall_at_45 =
      overall_true_in_top_45 /
      N_TRUE_POS,

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

    Overall_Coverage_80 =
      predictive_80$overall_coverage,

    Hidden_Coverage_80 =
      predictive_80$hidden_coverage,

    Overall_Interval_Score_80 =
      predictive_80$overall_interval_score,

    Hidden_Interval_Score_80 =
      predictive_80$hidden_interval_score,

    Overall_Coverage_95 =
      predictive_95$overall_coverage,

    Hidden_Coverage_95 =
      predictive_95$hidden_coverage,

    Overall_Interval_Score_95 =
      predictive_95$overall_interval_score,

    Hidden_Interval_Score_95 =
      predictive_95$hidden_interval_score,

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
      seq_len(N_NODES) %in%
        hidden_positive_idx
    ),

    positive_signal_group =
      POSITIVE_SIGNAL_GROUP,

    covariate_signal = as.integer(
      seq_len(N_NODES) %in%
        COVARIATE_SIGNAL_IDX
    ),

    network_signal = as.integer(
      seq_len(N_NODES) %in%
        NETWORK_SIGNAL_IDX
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

  node_scores$unlabeled_rank[
    unlabeled_idx
  ] <- unlabeled_rank

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
    positive_score_standardization = "within the 45 true positives",
    z_clip = SAR_Z_CLIP,
    log_hiding_weight = "-lambda * clipped standardized oracle covariate score",
    oracle_score_available_to_fitted_methods = FALSE
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
    predictive_interval_levels = PREDICTIVE_INTERVAL_LEVELS,
    n_predictive_draws = N_PREDICTIVE_DRAWS,
    interval_quantile_type = 1L,
    truth_for_coverage = "known latent TRUE_T; SAR-hidden positives count as T=1",
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

  ## All geometry and the SAR labels are regenerated from deterministic seeds
  ## when a partial checkpoint is resumed.
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

  ## ------------------------------------------------------------------------
  ## Replication-specific SAR hiding mechanism
  ## ------------------------------------------------------------------------

  sar_design <- build_sar_label_design(
    X_raw = covariate_set$X_raw,
    sar_seed = seed_row$sar_seed
  )

  observed_y <- sar_design$observed_y
  hidden_positive_idx <- sar_design$hidden_positive_idx
  observed_positive_idx <- sar_design$observed_positive_idx
  signal_type <- sar_design$signal_type

  sar_selection_summary_now <-
    sar_design$selection_summary

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
          c(
            "scenario",
            "replicate"
          )
        )
      ),
      drop = FALSE
    ]

  sar_selection_diagnostics_now <-
    sar_design$selection_diagnostics

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
          c(
            "scenario",
            "replicate"
          )
        )
      ),
      drop = FALSE
    ]

  cat(
    sprintf(
      paste0(
        "[%s] Replicate %03d/%03d | SAR hidden composition: ",
        "cov=%d, net=%d, both=%d | ",
        "mean score hidden=%.3f, observed=%.3f\n"
      ),
      timestamp_now(),
      replicate_number,
      N_REP,
      sar_selection_summary_now$hidden_covariates_only,
      sar_selection_summary_now$hidden_network_only,
      sar_selection_summary_now$hidden_both,
      sar_selection_summary_now$mean_score_hidden,
      sar_selection_summary_now$mean_score_observed
    )
  )

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
      sar_selection_summary = NULL,
      sar_selection_diagnostics = NULL,
      label_design = NULL,
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
      sar_seed = seed_row$sar_seed,
      lambda_SAR = LAMBDA_SAR,
      n = N_NODES,
      n_true_zero = N_TRUE_ZERO,
      n_true_positive = N_TRUE_POS,
      n_observed_positive = length(
        observed_positive_idx
      ),
      n_hidden_positive = length(
        hidden_positive_idx
      ),
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
  checkpoint$sar_selection_summary <-
    sar_selection_summary_now
  checkpoint$sar_selection_diagnostics <-
    sar_selection_diagnostics_now
  checkpoint$label_design <- list(
    observed_y = observed_y,
    hidden_positive_idx = hidden_positive_idx,
    observed_positive_idx = observed_positive_idx,
    signal_type = signal_type
  )

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
        Y = observed_y,
        hidden_positive_idx = hidden_positive_idx,
        observed_positive_idx = observed_positive_idx,
        positive_signal_group = POSITIVE_SIGNAL_GROUP,
        signal_type = signal_type,
        covariate_signal_idx = COVARIATE_SIGNAL_IDX,
        network_signal_idx = NETWORK_SIGNAL_IDX,
        covariate_profile =
          covariate_set$covariate_profile,
        X_raw = covariate_set$X_raw,
        X_std = X_cov,
        X_center = covariate_set$center,
        X_scale = covariate_set$scale,
        A_network = network_set$A,
        network_summary = network_set$summary,
        Z_network_modularity = Z_net,
        modularity_q = modularity$q,
        modularity_selected_eigenvalues =
          modularity$selected_eigenvalues,
        oracle_covariate_score =
          sar_design$oracle_score,
        oracle_covariate_block_scores =
          sar_design$oracle_block_scores,
        standardized_positive_oracle_score =
          sar_design$standardized_positive_score,
        positive_hidden_weights =
          sar_design$hidden_weights,
        sar_selection_summary =
          sar_selection_summary_now,
        sar_selection_diagnostics =
          sar_selection_diagnostics_now
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
        y = observed_y,
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

      saveRDS(
        checkpoint,
        checkpoint_file
      )

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
      replicate_number = replicate_number,
      observed_y = observed_y,
      hidden_positive_idx = hidden_positive_idx,
      signal_type = signal_type
    )

    evaluation$metrics$runtime_sec <-
      method_runtime_sec

    checkpoint$method_results[[method_id]] <-
      evaluation

    checkpoint$method_failures[[method_id]] <-
      NULL

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
    sar_design,
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

sar_selection_summary <- do.call(
  rbind,
  lapply(
    replication_results,
    `[[`,
    "sar_selection_summary"
  )
)

row.names(
  sar_selection_summary
) <- NULL

sar_selection_diagnostics <- do.call(
  rbind,
  lapply(
    replication_results,
    `[[`,
    "sar_selection_diagnostics"
  )
)

row.names(
  sar_selection_diagnostics
) <- NULL

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
    selection = paste0(
      "Exactly 15 of the 45 true positives are sampled without replacement ",
      "with weights proportional to exp(-lambda * z_i), where z_i is the ",
      "oracle covariate fraud-likeness score standardized among the true ",
      "positives and clipped to [-3,3]."
    ),
    oracle_score_used_for_fitting = FALSE
  ),

  sar_selection_summary = sar_selection_summary,
  sar_selection_Monte_Carlo_summary = sar_selection_mc_summary,

  methods = method_registry,

  metric_definitions = list(
    Overall_AUC = paste0(
      "AUC on all 400 observations against true T with 45 true positives."
    ),
    Hidden_AUC = paste0(
      "AUC among the 370 replication-specific Y=0 observations, containing ",
      "15 SAR-hidden positives and 355 true zeros."
    ),
    Overall_LogLoss =
      "LogLoss on all 400 observations against true T.",
    Hidden_LogLoss = paste0(
      "Mean -log(p_i) over the 15 replication-specific SAR-hidden positives, ",
      "all of which have T_i=1."
    ),
    Overall_Brier =
      "Brier score on all 400 observations against true T.",
    Hidden_Brier = paste0(
      "Mean (p_i-1)^2 over the 15 replication-specific SAR-hidden positives."
    ),
    Overall_Recall_at_45 = paste0(
      "Fraction of all 45 true positives whose descending midrank among ",
      "all 400 observations is at most 45."
    ),
    Hidden_Recall_at_15 = paste0(
      "Fraction of the 15 SAR-hidden positives whose descending midrank ",
      "among the 370 replication-specific Y=0 observations is at most 15."
    ),
    Coverage = paste0(
      "For each node, the central empirical binary predictive interval is ",
      "formed from predictive T draws using quantile type 1 and compared ",
      "with known TRUE_T. Overall averages over 400 nodes; Hidden averages ",
      "over the 15 replication-specific SAR-hidden positives."
    ),
    Interval_Score = paste0(
      "Equal-tail interval score evaluated against known TRUE_T, with ",
      "alpha=0.20 at 80% and alpha=0.05 at 95%."
    )
  ),

  ranking_rule =
    "rank(-probability, ties.method='average')",

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
    Ada =
      "NA: no model-based uncertainty in the implemented method.",
    Ada_Net =
      "NA: no model-based uncertainty in the implemented method.",
    GNN =
      "NA: no model-based uncertainty in the implemented method.",
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
  sar_selection_summary,
  sar_selection_diagnostics,
  sar_selection_mc_summary,
  failure_log,
  metric_summary_long,
  main_table_means,
  paper_table_mean_sd,
  uncertainty_table_mean_sd,
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
    "SIMULATION_100_COMPETING_METHODS_SAR_COMPLETE_RESULTS.RData"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    OUT_DIR,
    "sessionInfo_simulation_competing_methods_SAR.txt"
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
