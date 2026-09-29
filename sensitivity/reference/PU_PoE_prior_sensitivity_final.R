## ============================================================================
## PRIOR ROBUSTNESS ANALYSIS FOR THE GP-PU-POE TOY ILLUSTRATION
## ----------------------------------------------------------------------------
## This script reproduces exactly the seed-53 illustrative data and the
## GP-PU-PoE sampler used in the original four-model illustration. It then fits
## the same PU-PoE model under 11 distinct prior specifications:
##
##   1 default fit;
##   2 alternatives for the GP amplitudes;
##   2 alternatives for the GP lengthscales;
##   2 alternatives for the structured-mean scales;
##   2 alternatives for beta0;
##   2 alternatives for eta.
##
## Only one prior block is changed at a time. All fits use the same data, the
## same 10,000 / 20,000 / 150,000 / thin-30 MCMC schedule and the same MCMC
## random seed. The common seed is intentional: it holds the random-number
## stream fixed as far as possible while the prior specification changes.
##
## The final output is a 5 x 3 figure. The middle column repeats the one default
## fit in every row. The default colour is exactly green3, as in the original
## illustration; the conservative and diffuse columns use lighter and darker
## versions of the same colour.
## ============================================================================

options(stringsAsFactors = FALSE)

## ============================================================================
## 0. PACKAGE CHECKS
## ============================================================================

required_namespaces <- c(
  "Matrix",
  "igraph",
  "RSpectra",
  "coda",
  "ggplot2",
  "patchwork"
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

SCRIPT_VERSION <- "toy_seed53_PU_PoE_prior_sensitivity_HPD80_v1_0"

## Exact toy-data settings from the original illustration.
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

MODULARITY_SEED <- 53043L
MOD_EIG_TOL <- 1e-8
MOD_MAX_Q <- Inf
MOD_K_MAX <- 20L

## Exact MCMC schedule used by the default GP-PU-PoE illustration fit.
BLOCKED_BURNIN <- 10000L
JOINT_BURNIN <- 20000L
POSTERIOR_ITERATIONS <- 150000L
THIN <- 30L
N_RETAINED <- POSTERIOR_ITERATIONS %/% THIN
HPD_PROBABILITY <- 0.80

if (POSTERIOR_ITERATIONS %% THIN != 0L || N_RETAINED != 5000L) {
  stop("The requested sampling schedule must retain exactly 5,000 draws.")
}

## Baseline prior settings. The default fit uses these values exactly.
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

## MCMC proposal/adaptation controls copied from the original illustration.
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

## Ranking summaries saved with the figure.
TOP_K_VALUES <- c(10L, 20L, 30L)

## Exact original GP-PU-PoE MCMC seed. It is used for all 11 fits so that the
## prior specification is the only intended source of change.
COMMON_MCMC_SEED <- 530401L

## Plot colours.
BASE_PUPOE_COLOUR <- "green3"
CONSERVATIVE_LIGHTEN_AMOUNT <- 0.45
DIFFUSE_DARKEN_AMOUNT <- 0.28

OUT_DIR <- file.path(
  getwd(),
  "TOY_PU_POE_PRIOR_SENSITIVITY_SEED53"
)

FIT_CHECKPOINT_DIR <- file.path(
  OUT_DIR,
  "fit_checkpoints"
)

for (path_now in c(OUT_DIR, FIT_CHECKPOINT_DIR)) {
  dir.create(
    path_now,
    showWarnings = FALSE,
    recursive = TRUE
  )
}

RESUME_COMPLETED_FITS <- TRUE
FORCE_FRESH_RUN <- FALSE

if (isTRUE(FORCE_FRESH_RUN)) {
  unlink(
    FIT_CHECKPOINT_DIR,
    recursive = TRUE,
    force = TRUE
  )
  dir.create(
    FIT_CHECKPOINT_DIR,
    showWarnings = FALSE,
    recursive = TRUE
  )
}

RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

make_sensitivity_spec <- function(
    fit_id,
    varying_block,
    prior_role,
    total_sd_center = 2,
    corr_pct_range = c(0.10, 0.50),
    mean_rms_95_mult = 4.16,
    beta0_prev_range = c(0.01, 0.10),
    eta_prior_kind = "calibrated_q95",
    eta_alpha = 2,
    eta_beta = NA_real_,
    eta_q95 = 0.50,
    max_sigma_cov = 5,
    max_mean_sd_network = 5
) {
  list(
    fit_id = fit_id,
    varying_block = varying_block,
    prior_role = prior_role,
    total_sd_center = as.numeric(total_sd_center),
    corr_pct_range = as.numeric(corr_pct_range),
    mean_rms_95_mult = as.numeric(mean_rms_95_mult),
    beta0_prev_range = as.numeric(beta0_prev_range),
    eta_prior_kind = eta_prior_kind,
    eta_alpha = as.numeric(eta_alpha),
    eta_beta = as.numeric(eta_beta),
    eta_q95 = as.numeric(eta_q95),
    max_sigma_cov = as.numeric(max_sigma_cov),
    max_mean_sd_network = as.numeric(max_mean_sd_network),
    mcmc_seed = COMMON_MCMC_SEED
  )
}

## The conservative eta prior is Beta(1,5), which is the final choice agreed
## after the working-plan draft. The diffuse amplitude fit uses support 8 rather
## than 5 so that its requested central prior range up to 5.657 is not clipped.
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

SENSITIVITY_RUN_ORDER <- names(SENSITIVITY_SPECS)

BASE_PUPOE_MODEL_SPEC <- list(
  method = "GP-PU-PoE",
  active_experts = c("cov", "network"),
  use_pu = TRUE,
  mcmc_seed = COMMON_MCMC_SEED
)

## ============================================================================
## 2. GENERAL HELPERS
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
  list(
    scaled = scaled,
    center = center,
    scale = scale_value
  )
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
    if (!is.null(R)) {
      return(list(R = R, jitter = jj))
    }
  }
  stop("Cholesky factorization failed after jitter escalation.")
}

svd_psd_fallback <- function(K) {
  K <- symmetrize(K)
  E <- tryCatch(
    eigen(K, symmetric = TRUE),
    error = function(e) NULL
  )
  if (!is.null(E)) {
    return(list(
      U = Re(E$vectors),
      d = pmax(Re(E$values), 0)
    ))
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

rInvGamma <- function(n, shape, scale) {
  1 / rgamma(n, shape = shape, rate = scale)
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

round_prior <- function(x) {
  round(as.numeric(x), PRIOR_ROUND_DIGITS)
}

logit <- function(p) {
  log(p / (1 - p))
}

file_md5_or_na <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path)[1L])
}

## ============================================================================
## 3. GENERATE THE EXACT USER-SUPPLIED TOY DATA
## ============================================================================

generate_exact_toy_seed53 <- function() {
  ## The statements below intentionally mirror the supplied generator.
  N_NODES_LOCAL <- 100L
  N_COV_LOCAL <- 10L
  POS_IDX_LOCAL <- 71:100
  OBSERVED_IDX_LOCAL <- 81:100
  HIDDEN_IDX_LOCAL <- 71:80
  ZERO_IDX_LOCAL <- 1:70
  SIGNAL_MEAN_LOCAL <- 0.5
  BG_MEAN_LOCAL <- 0.0
  COMMON_SD_LOCAL <- 1.0
  P00_LOCAL <- 0.10
  P01_LOCAL <- 0.05
  P11_LOCAL <- 0.16

  set.seed(53)

  X <- matrix(
    rnorm(
      N_NODES_LOCAL * N_COV_LOCAL,
      BG_MEAN_LOCAL,
      COMMON_SD_LOCAL
    ),
    N_NODES_LOCAL,
    N_COV_LOCAL
  )

  X[POS_IDX_LOCAL, ] <- rnorm(
    length(POS_IDX_LOCAL) * N_COV_LOCAL,
    SIGNAL_MEAN_LOCAL,
    COMMON_SD_LOCAL
  )

  colnames(X) <- paste0("X", seq_len(N_COV_LOCAL))

  Pmat <- matrix(P00_LOCAL, N_NODES_LOCAL, N_NODES_LOCAL)
  Pmat[POS_IDX_LOCAL, ] <- P01_LOCAL
  Pmat[, POS_IDX_LOCAL] <- P01_LOCAL
  Pmat[POS_IDX_LOCAL, POS_IDX_LOCAL] <- P11_LOCAL

  A <- matrix(0L, N_NODES_LOCAL, N_NODES_LOCAL)
  ut <- upper.tri(A)
  A[ut] <- as.integer(runif(sum(ut)) < Pmat[ut])
  A <- A + t(A)
  diag(A) <- 0L

  T <- integer(N_NODES_LOCAL)
  T[POS_IDX_LOCAL] <- 1L

  Y <- integer(N_NODES_LOCAL)
  Y[OBSERVED_IDX_LOCAL] <- 1L

  list(
    seed = 53L,
    N = N_NODES_LOCAL,
    X = X,
    A = A,
    T = T,
    Y = Y,
    pos_idx = POS_IDX_LOCAL,
    observed_idx = OBSERVED_IDX_LOCAL,
    hidden_idx = HIDDEN_IDX_LOCAL,
    zero_idx = ZERO_IDX_LOCAL,
    probability_matrix = Pmat
  )
}

toy <- generate_exact_toy_seed53()

if (
  toy$seed != TOY_SEED ||
  toy$N != N_NODES ||
  !identical(as.integer(toy$zero_idx), TRUE_ZERO_IDX) ||
  !identical(as.integer(toy$hidden_idx), HIDDEN_POS_IDX) ||
  !identical(as.integer(toy$observed_idx), OBSERVED_POS_IDX) ||
  !identical(as.integer(toy$pos_idx), TRUE_POS_IDX) ||
  !identical(which(toy$T == 1L), TRUE_POS_IDX) ||
  !identical(which(toy$Y == 1L), OBSERVED_POS_IDX) ||
  !identical(which(toy$T == 1L & toy$Y == 0L), HIDDEN_POS_IDX)
) {
  stop("The exact toy label design was not reproduced.")
}

if (
  !all(dim(toy$X) == c(N_NODES, N_COVARIATES)) ||
  !all(dim(toy$A) == c(N_NODES, N_NODES)) ||
  any(!is.finite(toy$X)) ||
  any(toy$A != t(toy$A)) ||
  any(diag(toy$A) != 0L) ||
  any(!toy$A %in% c(0L, 1L))
) {
  stop("The exact toy X or A object is invalid.")
}

TOY_RDS_FILE <- file.path(OUT_DIR, "toy_seed53.rds")
TOY_COVARIATE_FILE <- file.path(OUT_DIR, "toy_seed53_covariates.csv")
TOY_ADJACENCY_FILE <- file.path(OUT_DIR, "toy_seed53_adjacency.csv")

saveRDS(toy, TOY_RDS_FILE)

write.csv(
  cbind(
    node = seq_len(N_NODES),
    T = toy$T,
    Y = toy$Y,
    toy$X
  ),
  TOY_COVARIATE_FILE,
  row.names = FALSE
)

write.csv(
  toy$A,
  TOY_ADJACENCY_FILE,
  row.names = FALSE
)

TOY_DATA_HASH <- file_md5_or_na(TOY_RDS_FILE)

## Model covariates are centered and sample-standardized, as in the supplied
## simulation/application model code. The raw generated X is retained above.
X_standardization <- standardize_columns(toy$X)
X_cov <- X_standardization$scaled
colnames(X_cov) <- colnames(toy$X)
A_network <- sanitize_network_adjacency(toy$A)
Y <- as.integer(toy$Y)
T_TRUE <- as.integer(toy$T)

NODE_GROUP <- rep("true zero", N_NODES)
NODE_GROUP[HIDDEN_POS_IDX] <- "hidden positive"
NODE_GROUP[OBSERVED_POS_IDX] <- "observed positive"

## ============================================================================
## 4. LABEL-BLIND MODULARITY FEATURES
## ============================================================================

choose_q_modularity_eigengap <- function(
    values,
    tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q
) {
  positive_values <- sort(
    as.numeric(values[values > tol]),
    decreasing = TRUE
  )

  n_positive <- length(positive_values)

  if (n_positive == 0L) {
    return(list(
      q = 1L,
      positive_values = numeric(0),
      gaps = numeric(0),
      max_gap = NA_real_
    ))
  }

  if (n_positive == 1L) {
    return(list(
      q = 1L,
      positive_values = positive_values,
      gaps = numeric(0),
      max_gap = NA_real_
    ))
  }

  gaps <- positive_values[-n_positive] - positive_values[-1L]
  q <- which.max(gaps)

  if (is.finite(max_q)) {
    q <- min(q, as.integer(max_q))
  }

  q <- max(1L, min(q, n_positive))

  list(
    q = as.integer(q),
    positive_values = positive_values,
    gaps = gaps,
    max_gap = gaps[which.max(gaps)]
  )
}

orient_eigenvectors_deterministically <- function(V) {
  V <- as.matrix(V)
  for (column_index in seq_len(ncol(V))) {
    largest_loading_index <- which.max(abs(V[, column_index]))
    if (V[largest_loading_index, column_index] < 0) {
      V[, column_index] <- -V[, column_index]
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
  A <- sanitize_network_adjacency(A)
  n <- nrow(A)
  degree <- rowSums(A)
  n_edges <- sum(degree) / 2

  if (!is.finite(n_edges) || n_edges <= 0) {
    stop("Cannot construct modularity features from an edgeless network.")
  }

  modularity_operator <- function(x, args) {
    as.vector(
      A %*% x -
        degree * sum(degree * x) / (2 * n_edges)
    )
  }

  k_use <- min(as.integer(k_max), n - 2L)

  eig <- RSpectra::eigs_sym(
    A = modularity_operator,
    k = k_use,
    n = n,
    which = "LA",
    opts = list(retvec = TRUE)
  )

  eigenvalues <- Re(eig$values)
  eigenvectors <- Re(eig$vectors)
  ordering <- order(eigenvalues, decreasing = TRUE)
  eigenvalues <- eigenvalues[ordering]
  eigenvectors <- eigenvectors[, ordering, drop = FALSE]
  eigenvectors <- orient_eigenvectors_deterministically(eigenvectors)

  positive_index <- which(eigenvalues > eig_tol)

  if (length(positive_index) == 0L) {
    keep <- 1L
  } else {
    gap_information <- choose_q_modularity_eigengap(
      values = eigenvalues,
      tol = eig_tol,
      max_q = max_q
    )
    keep <- positive_index[seq_len(gap_information$q)]
  }

  Z_raw <- eigenvectors[, keep, drop = FALSE]
  colnames(Z_raw) <- paste0("network_mod", seq_len(ncol(Z_raw)))

  standardized <- standardize_columns(Z_raw)
  Z <- standardized$scaled
  colnames(Z) <- colnames(Z_raw)

  list(
    Z = Z,
    raw = Z_raw,
    q = ncol(Z),
    selected_eigenvalues = eigenvalues[keep],
    all_eigenvalues = eigenvalues,
    center = standardized$center,
    scale = standardized$scale
  )
}

modularity <- modularity_features_single_network(
  A = A_network,
  seed = MODULARITY_SEED,
  eig_tol = MOD_EIG_TOL,
  max_q = MOD_MAX_Q,
  k_max = MOD_K_MAX
)

Z_network <- as.matrix(modularity$Z)

if (
  nrow(Z_network) != N_NODES ||
  ncol(Z_network) < 1L ||
  any(!is.finite(Z_network))
) {
  stop("The toy modularity feature matrix is invalid.")
}

write.csv(
  data.frame(
    node = seq_len(N_NODES),
    Z_network,
    check.names = FALSE
  ),
  file.path(OUT_DIR, "toy_seed53_modularity_features.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    q_network = modularity$q,
    selected_eigenvalue = modularity$selected_eigenvalues
  ),
  file.path(OUT_DIR, "toy_seed53_modularity_eigenvalues.csv"),
  row.names = FALSE
)

cat("\n================ EXACT TOY DATA CHECK ================\n")
cat("Seed                    :", TOY_SEED, "\n")
cat("Nodes                   :", N_NODES, "\n")
cat("True zeros              :", length(TRUE_ZERO_IDX), "\n")
cat("Hidden positives        :", length(HIDDEN_POS_IDX), "\n")
cat("Observed positives      :", length(OBSERVED_POS_IDX), "\n")
cat("Covariates              :", ncol(X_cov), "\n")
cat("Network edges           :", sum(A_network) / 2, "\n")
cat("Modularity features q   :", ncol(Z_network), "\n")
cat("SBM p00/p01/p11         :", P_00, "/", P_01, "/", P_11, "\n")
cat("Toy-data MD5            :", TOY_DATA_HASH, "\n")


## ============================================================================
## 5. EXACT NETWORK-GP GEOMETRY
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

network_geometry <- prepare_network_components_for_blocks(
  A_network,
  network_name = "toy_network"
)

write.csv(
  network_geometry$component_summary,
  file.path(OUT_DIR, "toy_seed53_network_component_summary.csv"),
  row.names = FALSE
)

## ============================================================================
## 6. FIXED BETA0 AND ETA PRIORS
## ============================================================================

solve_eta_beta_from_q95 <- function(
    alpha,
    q95_target,
    prob = 0.95
) {
  if (
    !is.finite(alpha) || alpha <= 0 ||
    !is.finite(q95_target) ||
    q95_target <= 0 ||
    q95_target >= 1
  ) {
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

if (
  !is.finite(beta0_prior_for_sampler$mean) ||
  !is.finite(beta0_prior_for_sampler$sd) ||
  beta0_prior_for_sampler$sd <= 0
) {
  stop("The beta0 prior is invalid.")
}

if (
  !is.finite(eta_prior_for_sampler$a_eta) ||
  !is.finite(eta_prior_for_sampler$b_eta) ||
  eta_prior_for_sampler$a_eta <= 0 ||
  eta_prior_for_sampler$b_eta <= 0
) {
  stop("The eta prior is invalid.")
}

## ============================================================================
## 7. GEOMETRY-SPECIFIC PRIOR CALIBRATION
## ============================================================================

fit_lognormal_interval <- function(
    q_low,
    q_high,
    probs = c(0.025, 0.975)
) {
  q_low <- as.numeric(q_low)
  q_high <- as.numeric(q_high)

  if (
    length(q_low) != 1L ||
    length(q_high) != 1L ||
    !is.finite(q_low) ||
    !is.finite(q_high) ||
    q_low <= 0 ||
    q_high <= 0
  ) {
    stop("q_low and q_high must be positive finite scalars.")
  }

  qlo <- min(q_low, q_high)
  qhi <- max(q_low, q_high)
  p_low <- probs[1L]
  p_high <- probs[2L]

  if (
    !is.finite(p_low) ||
    !is.finite(p_high) ||
    p_low <= 0 ||
    p_high >= 1 ||
    p_low >= p_high
  ) {
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

    if (
      length(ell) != 1L ||
      !is.finite(ell) ||
      ell <= 0
    ) {
      return(list(
        mean_corr = NA_real_,
        mean_sd = NA_real_
      ))
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

    if (
      length(ell) != 1L ||
      !is.finite(ell) ||
      ell <= 0
    ) {
      return(list(
        mean_corr = NA_real_,
        mean_sd = NA_real_
      ))
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

      if (
        is.finite(off_corr_sum) &&
        off_corr_sum < 0 &&
        abs(off_corr_sum) < 1e-8
      ) {
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
    is.finite(curve$ell) &
      is.finite(curve$mean_corr),
    ,
    drop = FALSE
  ]

  curve <- curve[order(curve$ell), , drop = FALSE]

  corr_range <- range(curve$mean_corr, na.rm = TRUE)
  target_used <- min(
    max(target_corr, corr_range[1L]),
    corr_range[2L]
  )

  difference <- curve$mean_corr - target_used
  exact_idx <- which(abs(difference) <= 1e-12)

  if (length(exact_idx) > 0L) {
    ell_hat <- curve$ell[exact_idx[1L]]
  } else {
    crossing_idx <- which(
      difference[-nrow(curve)] *
        difference[-1L] <= 0
    )

    root_result <- NULL

    if (length(crossing_idx) > 0L) {
      crossing_score <- abs(difference[crossing_idx]) +
        abs(difference[crossing_idx + 1L])
      kk <- crossing_idx[which.min(crossing_score)]
      log_interval <- log(c(
        curve$ell[kk],
        curve$ell[kk + 1L]
      ))

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

calibrate_one_expert_prior <- function(
    expert_name,
    metric_fun,
    ell_bounds,
    effective_sd_center,
    corr_pct_range = CORR_PCT_RANGE_ALL,
    sigma_range_factor = SIGMA_RANGE_FACTOR
) {
  corr_pct_range <- as.numeric(corr_pct_range)
  sigma_range_factor <- as.numeric(sigma_range_factor)

  if (
    length(corr_pct_range) != 2L ||
    any(!is.finite(corr_pct_range)) ||
    corr_pct_range[1L] < 0 ||
    corr_pct_range[2L] > 1 ||
    corr_pct_range[1L] >= corr_pct_range[2L]
  ) {
    stop("corr_pct_range must be an increasing pair inside [0,1].")
  }

  if (
    length(sigma_range_factor) != 1L ||
    !is.finite(sigma_range_factor) ||
    sigma_range_factor <= 1
  ) {
    stop("sigma_range_factor must be one finite value larger than one.")
  }

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
    corr_pct_range * attainable_span

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
    effective_sd_center / sigma_range_factor
  ) / unit_kernel_mean_sd

  sigma_raw_high <- (
    effective_sd_center * sigma_range_factor
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
    target_corr_pct_low = corr_pct_range[1L],
    target_corr_pct_high = corr_pct_range[2L],
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
    effective_sd_target_low = effective_sd_center / sigma_range_factor,
    effective_sd_target_high = effective_sd_center * sigma_range_factor,
    sigma_range_factor = sigma_range_factor,
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

  if (output$log_ell_sd <= 0) {
    output$log_ell_sd <- hyperprior_unrounded$log_ell_sd
  }

  if (output$log_sig_sd <= 0) {
    output$log_sig_sd <- hyperprior_unrounded$log_sig_sd
  }

  if (
    any(!is.finite(unlist(output))) ||
    output$log_ell_sd <= 0 ||
    output$log_sig_sd <= 0
  ) {
    stop("The active geometry-specific hyperprior is invalid.")
  }

  output
}

calibrate_model_priors <- function(
    model_spec,
    X_cov,
    Z_network,
    network_geometry
) {
  active_experts <- model_spec$active_experts
  n_experts <- length(active_experts)

  if (n_experts < 1L || n_experts > 2L) {
    stop("The toy illustration supports one or two active experts.")
  }

  ## A one-expert ablation gets the full total SD budget. In a two-expert model
  ## the two experts split total residual variance equally, exactly as in the
  ## supplied full model.
  effective_sd_center <- if (n_experts == 1L) {
    setNames(TOTAL_SD_CENTER, active_experts)
  } else {
    setNames(
      rep(TOTAL_SD_CENTER / sqrt(n_experts), n_experts),
      active_experts
    )
  }

  active_hyperpriors <- list()
  hyperpriors_unrounded <- list()
  summary_rows <- list()
  curves <- list()
  tau_scale_unrounded <- numeric(n_experts)
  names(tau_scale_unrounded) <- active_experts

  half_t_q95_unit <- qt(
    (0.95 + 1) / 2,
    df = HALFT_DF
  )

  for (expert_name in active_experts) {
    if (identical(expert_name, "cov")) {
      metric_fun <- make_covariate_metric_fun(X_cov)
      ell_bounds <- ELL_BOUNDS_COV
      mean_design <- X_cov
    } else if (identical(expert_name, "network")) {
      metric_fun <- network_metric_fun_factory(network_geometry)
      ell_bounds <- ELL_BOUNDS_NETWORK
      mean_design <- Z_network
    } else {
      stop("Unknown expert: ", expert_name)
    }

    calibration <- calibrate_one_expert_prior(
      expert_name = expert_name,
      metric_fun = metric_fun,
      ell_bounds = ell_bounds,
      effective_sd_center = effective_sd_center[expert_name]
    )

    active_hyperpriors[[expert_name]] <- activate_hyperprior(
      calibration$hyperprior_unrounded
    )

    hyperpriors_unrounded[[expert_name]] <-
      calibration$hyperprior_unrounded

    b_factor <- sqrt(
      sum(mean_design^2) /
        nrow(mean_design)
    )

    if (!is.finite(b_factor) || b_factor <= 0) {
      stop("The structured-mean RMS factor is invalid for ", expert_name, ".")
    }

    tau_scale_unrounded[expert_name] <- (
      MEAN_RMS_95_MULT * effective_sd_center[expert_name]
    ) / (
      half_t_q95_unit * b_factor
    )

    row <- calibration$row
    row$model_id <- model_spec$model_id
    row$method <- model_spec$method
    row$n_active_experts <- n_experts
    row$n_mean_features <- ncol(mean_design)
    row$Z_rms_factor <- b_factor
    row$tau_scale_unrounded <- tau_scale_unrounded[expert_name]
    summary_rows[[expert_name]] <- row
    curves[[expert_name]] <- calibration$curve
  }

  tau_scale <- if (isTRUE(ROUND_ACTIVE_PRIORS)) {
    round_prior(tau_scale_unrounded)
  } else {
    tau_scale_unrounded
  }
  names(tau_scale) <- active_experts

  if (any(!is.finite(tau_scale)) || any(tau_scale <= 0)) {
    stop("At least one active half-t scale is invalid.")
  }

  tau2_default <- tau_scale^2

  summary_table <- do.call(rbind, summary_rows)
  row.names(summary_table) <- NULL
  summary_table$tau_scale_active <- unname(
    tau_scale[summary_table$expert]
  )
  summary_table$tau2_default <- unname(
    tau2_default[summary_table$expert]
  )
  summary_table$log_ell_mean_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_ell_mean
    },
    numeric(1)
  )
  summary_table$log_ell_sd_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_ell_sd
    },
    numeric(1)
  )
  summary_table$log_sig_mean_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_sig_mean
    },
    numeric(1)
  )
  summary_table$log_sig_sd_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_sig_sd
    },
    numeric(1)
  )
  summary_table$active_prior_rounding <- ROUND_ACTIVE_PRIORS
  summary_table$active_prior_digits <- PRIOR_ROUND_DIGITS

  list(
    model_id = model_spec$model_id,
    method = model_spec$method,
    active_experts = active_experts,
    effective_sd_center = effective_sd_center,
    active_hyperpriors = active_hyperpriors,
    hyperpriors_unrounded = hyperpriors_unrounded,
    tau_scale = tau_scale,
    tau_scale_unrounded = tau_scale_unrounded,
    tau2_default = tau2_default,
    summary = summary_table,
    curves = curves
  )
}


## ============================================================================
## 7A. BUILD THE 11 PRIOR BUNDLES
## ============================================================================

make_beta0_prior_from_range <- function(prevalence_range) {
  prevalence_range <- as.numeric(prevalence_range)

  if (
    length(prevalence_range) != 2L ||
    any(!is.finite(prevalence_range)) ||
    prevalence_range[1L] <= 0 ||
    prevalence_range[2L] >= 1 ||
    prevalence_range[1L] >= prevalence_range[2L]
  ) {
    stop("The beta0 prevalence range is invalid.")
  }

  lower_logit <- logit(prevalence_range[1L])
  upper_logit <- logit(prevalence_range[2L])

  list(
    mean = round_prior(0.5 * (lower_logit + upper_logit)),
    sd = round_prior(
      (upper_logit - lower_logit) /
        (2 * qnorm(0.975))
    ),
    prevalence_range = prevalence_range
  )
}

make_eta_prior_from_spec <- function(fit_spec) {
  if (identical(fit_spec$eta_prior_kind, "calibrated_q95")) {
    beta_value <- solve_eta_beta_from_q95(
      alpha = fit_spec$eta_alpha,
      q95_target = fit_spec$eta_q95,
      prob = 0.95
    )

    return(list(
      a_eta = as.numeric(fit_spec$eta_alpha),
      b_eta = round_prior(beta_value),
      construction = "calibrated_q95",
      q95_target = as.numeric(fit_spec$eta_q95)
    ))
  }

  if (identical(fit_spec$eta_prior_kind, "fixed_beta")) {
    if (
      !is.finite(fit_spec$eta_alpha) ||
      !is.finite(fit_spec$eta_beta) ||
      fit_spec$eta_alpha <= 0 ||
      fit_spec$eta_beta <= 0
    ) {
      stop("The fixed Beta prior for eta is invalid.")
    }

    return(list(
      a_eta = as.numeric(fit_spec$eta_alpha),
      b_eta = as.numeric(fit_spec$eta_beta),
      construction = "fixed_beta",
      q95_target = NA_real_
    ))
  }

  stop("Unknown eta_prior_kind: ", fit_spec$eta_prior_kind)
}

calibrate_sensitivity_fit_priors <- function(
    fit_spec,
    X_cov,
    Z_network,
    network_geometry
) {
  active_experts <- c("cov", "network")
  n_experts <- length(active_experts)

  effective_sd_center <- setNames(
    rep(fit_spec$total_sd_center / sqrt(n_experts), n_experts),
    active_experts
  )

  ## The structured-mean scale is a separate prior block. It is therefore
  ## referenced to the default expert SD budget even in the sigma-sensitivity
  ## fits, so changing TOTAL_SD_CENTER changes only the GP amplitude priors.
  structured_mean_sd_reference <- setNames(
    rep(TOTAL_SD_CENTER / sqrt(n_experts), n_experts),
    active_experts
  )

  active_hyperpriors <- list()
  hyperpriors_unrounded <- list()
  summary_rows <- list()
  curves <- list()
  tau_scale_unrounded <- numeric(n_experts)
  names(tau_scale_unrounded) <- active_experts

  half_t_q95_unit <- qt(
    (0.95 + 1) / 2,
    df = HALFT_DF
  )

  for (expert_name in active_experts) {
    if (identical(expert_name, "cov")) {
      metric_fun <- make_covariate_metric_fun(X_cov)
      ell_bounds <- ELL_BOUNDS_COV
      mean_design <- X_cov
    } else {
      metric_fun <- network_metric_fun_factory(network_geometry)
      ell_bounds <- ELL_BOUNDS_NETWORK
      mean_design <- Z_network
    }

    calibration <- calibrate_one_expert_prior(
      expert_name = expert_name,
      metric_fun = metric_fun,
      ell_bounds = ell_bounds,
      effective_sd_center = effective_sd_center[expert_name],
      corr_pct_range = fit_spec$corr_pct_range,
      sigma_range_factor = SIGMA_RANGE_FACTOR
    )

    active_hyperpriors[[expert_name]] <- activate_hyperprior(
      calibration$hyperprior_unrounded
    )

    hyperpriors_unrounded[[expert_name]] <-
      calibration$hyperprior_unrounded

    b_factor <- sqrt(
      sum(mean_design^2) /
        nrow(mean_design)
    )

    if (!is.finite(b_factor) || b_factor <= 0) {
      stop("The structured-mean RMS factor is invalid for ", expert_name, ".")
    }

    tau_scale_unrounded[expert_name] <- (
      fit_spec$mean_rms_95_mult * structured_mean_sd_reference[expert_name]
    ) / (
      half_t_q95_unit * b_factor
    )

    row <- calibration$row
    row$fit_id <- fit_spec$fit_id
    row$varying_block <- fit_spec$varying_block
    row$prior_role <- fit_spec$prior_role
    row$n_active_experts <- n_experts
    row$n_mean_features <- ncol(mean_design)
    row$Z_rms_factor <- b_factor
    row$mean_rms_95_mult <- fit_spec$mean_rms_95_mult
    row$total_sd_center <- fit_spec$total_sd_center
    row$structured_mean_sd_reference <-
      structured_mean_sd_reference[expert_name]
    row$tau_scale_unrounded <- tau_scale_unrounded[expert_name]
    summary_rows[[expert_name]] <- row
    curves[[expert_name]] <- calibration$curve
  }

  tau_scale <- if (isTRUE(ROUND_ACTIVE_PRIORS)) {
    round_prior(tau_scale_unrounded)
  } else {
    tau_scale_unrounded
  }
  names(tau_scale) <- active_experts

  if (any(!is.finite(tau_scale)) || any(tau_scale <= 0)) {
    stop("At least one active half-t scale is invalid.")
  }

  tau2_default <- tau_scale^2
  beta0_prior <- make_beta0_prior_from_range(
    fit_spec$beta0_prev_range
  )
  eta_prior <- make_eta_prior_from_spec(fit_spec)

  summary_table <- do.call(rbind, summary_rows)
  row.names(summary_table) <- NULL

  summary_table$tau_scale_active <- unname(
    tau_scale[summary_table$expert]
  )
  summary_table$tau2_default <- unname(
    tau2_default[summary_table$expert]
  )
  summary_table$log_ell_mean_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_ell_mean
    },
    numeric(1)
  )
  summary_table$log_ell_sd_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_ell_sd
    },
    numeric(1)
  )
  summary_table$log_sig_mean_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_sig_mean
    },
    numeric(1)
  )
  summary_table$log_sig_sd_active <- vapply(
    summary_table$expert,
    function(expert_name) {
      active_hyperpriors[[expert_name]]$log_sig_sd
    },
    numeric(1)
  )
  summary_table$beta0_prior_mean <- beta0_prior$mean
  summary_table$beta0_prior_sd <- beta0_prior$sd
  summary_table$beta0_prev_low <- beta0_prior$prevalence_range[1L]
  summary_table$beta0_prev_high <- beta0_prior$prevalence_range[2L]
  summary_table$eta_prior_alpha <- eta_prior$a_eta
  summary_table$eta_prior_beta <- eta_prior$b_eta
  summary_table$max_sigma_cov <- fit_spec$max_sigma_cov
  summary_table$max_mean_sd_network <- fit_spec$max_mean_sd_network
  summary_table$active_prior_rounding <- ROUND_ACTIVE_PRIORS
  summary_table$active_prior_digits <- PRIOR_ROUND_DIGITS

  list(
    fit_id = fit_spec$fit_id,
    fit_spec = fit_spec,
    active_experts = active_experts,
    effective_sd_center = effective_sd_center,
    structured_mean_sd_reference = structured_mean_sd_reference,
    active_hyperpriors = active_hyperpriors,
    hyperpriors_unrounded = hyperpriors_unrounded,
    tau_scale = tau_scale,
    tau_scale_unrounded = tau_scale_unrounded,
    tau2_default = tau2_default,
    beta0_prior = beta0_prior,
    eta_prior = eta_prior,
    summary = summary_table,
    curves = curves
  )
}

sensitivity_prior_bundles <- lapply(
  SENSITIVITY_SPECS,
  calibrate_sensitivity_fit_priors,
  X_cov = X_cov,
  Z_network = Z_network,
  network_geometry = network_geometry
)

names(sensitivity_prior_bundles) <- SENSITIVITY_RUN_ORDER

prior_calibration_table <- do.call(
  rbind,
  lapply(sensitivity_prior_bundles, `[[`, "summary")
)
row.names(prior_calibration_table) <- NULL

write.csv(
  prior_calibration_table,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_calibration.csv"
  ),
  row.names = FALSE
)

cat("\n================ PRIOR SENSITIVITY FITS ================\n")
cat("Distinct fits          :", length(SENSITIVITY_RUN_ORDER), "\n")
cat("Common MCMC seed       :", COMMON_MCMC_SEED, "\n")
cat("Default beta0 prior    : N(",
    sensitivity_prior_bundles$DEFAULT$beta0_prior$mean, ", ",
    sensitivity_prior_bundles$DEFAULT$beta0_prior$sd, "^2)\n", sep = "")
cat("Default eta prior      : Beta(",
    sensitivity_prior_bundles$DEFAULT$eta_prior$a_eta, ", ",
    sensitivity_prior_bundles$DEFAULT$eta_prior$b_eta, ")\n", sep = "")
cat("Conservative eta prior : Beta(1, 5)\n")
cat("80% node intervals     : coda::HPDinterval\n")

## ============================================================================
## 8. PROPOSAL-COVARIANCE HELPERS
## ============================================================================

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

proposal_cov_scaled_2x2 <- function(
    kappa_scalar,
    prior_sd,
    cor_mat,
    jitter = 1e-12
) {
  kappa_scalar <- as.numeric(kappa_scalar)
  prior_sd <- as.numeric(prior_sd)

  if (
    length(kappa_scalar) != 1L ||
    !is.finite(kappa_scalar) ||
    kappa_scalar <= 0
  ) {
    stop("kappa_scalar must be one positive finite value.")
  }

  if (
    length(prior_sd) != 2L ||
    any(!is.finite(prior_sd)) ||
    any(prior_sd <= 0)
  ) {
    stop("prior_sd must contain two positive finite values.")
  }

  R <- sanitize_cor_2x2(cor_mat, max_abs = 0.999)
  D <- diag(prior_sd, 2L)
  S <- kappa_scalar * D %*% R %*% D
  symmetrize(S) + diag(jitter, 2L)
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

## ============================================================================
## 9. BERNOULLI AND PU LIKELIHOODS
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
    dbeta(eta, a_eta, b_eta, log = TRUE) +
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
  u0 <- qlogis(pmin(pmax(eta, eps), 1 - eps))

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
    if (!is.finite(value_left) || value_left <= log_height) break
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
    if (!is.finite(value_right) || value_right <= log_height) break
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

    if (is.finite(value_1) && value_1 >= log_height) break
    if (u1 < u0) left <- u1 else right <- u1
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

hidden_probability_given_zero <- function(p_latent, eta) {
  p_latent <- pmin(pmax(as.numeric(p_latent), 0), 1)
  eta <- pmin(pmax(as.numeric(eta), ETA_EPS), 1 - ETA_EPS)
  denominator <- 1 - (1 - eta) * p_latent
  denominator <- pmax(denominator, .Machine$double.eps)
  output <- eta * p_latent / denominator
  pmin(pmax(output, 0), 1)
}

## ============================================================================
## 10. GENERIC EXACT GP SAMPLER FOR ONE OR TWO ACTIVE EXPERTS
## ============================================================================

run_exact_toy_model_chain <- function(
    model_spec,
    prior_bundle,
    Y,
    X_cov,
    A_network,
    Z_network,
    network_geometry,
    m_beta,
    s_beta,
    a_eta,
    b_eta,
    n_blocked_burnin = BLOCKED_BURNIN,
    n_joint_burnin = JOINT_BURNIN,
    n_sample = POSTERIOR_ITERATIONS,
    thin = THIN,
    delta_init = DELTA_INIT,
    kappa_init = KAPPA_INIT,
    adapt_block_cor = ADAPT_BLOCK_CORRELATIONS,
    cor_adapt_start = COR_ADAPT_START,
    cor_adapt_interval = COR_ADAPT_INTERVAL,
    cor_adapt_window = COR_ADAPT_WINDOW,
    cor_shrinkage = COR_SHRINKAGE,
    cor_max_abs = COR_MAX_ABS,
    joint_cor_shrinkage = JOINT_COR_SHRINKAGE,
    joint_cor_max_abs = JOINT_COR_MAX_ABS,
    f_refresh_after_joint = F_REFRESH_AFTER_JOINT,
    eta_slice_width = ETA_SLICE_WIDTH,
    initial_eta = INITIAL_ETA,
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
    blocked_progress_every = BLOCKED_PROGRESS_EVERY,
    joint_progress_every = JOINT_PROGRESS_EVERY,
    sample_progress_every = SAMPLE_PROGRESS_EVERY,
    verbose = TRUE
) {
  Y <- as.numeric(Y)
  X_cov <- as.matrix(X_cov)
  A_network <- sanitize_network_adjacency(A_network)
  Z_network <- as.matrix(Z_network)

  active_experts <- as.character(model_spec$active_experts)
  use_pu <- isTRUE(model_spec$use_pu)
  n_experts <- length(active_experts)
  n <- length(Y)

  n_blocked_burnin <- as.integer(n_blocked_burnin)
  n_joint_burnin <- as.integer(n_joint_burnin)
  n_sample <- as.integer(n_sample)
  thin <- as.integer(thin)
  adapt_start <- as.integer(adapt_start)
  adapt_interval <- as.integer(adapt_interval)
  cor_adapt_start <- as.integer(cor_adapt_start)
  cor_adapt_interval <- as.integer(cor_adapt_interval)
  cor_adapt_window <- as.integer(cor_adapt_window)
  emp_cov_tail <- as.integer(emp_cov_tail)

  if (!all(Y %in% c(0, 1))) {
    stop("Y must contain only 0 and 1.")
  }

  if (
    nrow(X_cov) != n ||
    nrow(Z_network) != n ||
    !all(dim(A_network) == c(n, n))
  ) {
    stop("The model inputs have incompatible dimensions.")
  }

  if (
    n_blocked_burnin != 10000L ||
    n_joint_burnin != 20000L ||
    n_sample != 150000L ||
    thin != 30L
  ) {
    stop("The requested 10k/20k/150k/thin-30 schedule was altered.")
  }

  if (n_sample %% thin != 0L) {
    stop("n_sample must be divisible by thin.")
  }

  if (
    n_experts < 1L ||
    n_experts > 2L ||
    any(!active_experts %in% c("cov", "network")) ||
    anyDuplicated(active_experts)
  ) {
    stop("active_experts must be a unique subset of cov and network.")
  }

  ## ------------------------------------------------------------------------
  ## Mean design and expert-index bookkeeping
  ## ------------------------------------------------------------------------

  mean_design_list <- list()

  for (expert_name in active_experts) {
    if (identical(expert_name, "cov")) {
      mean_design_list[[expert_name]] <- X_cov
    } else {
      mean_design_list[[expert_name]] <- Z_network
    }
  }

  W_mean <- do.call(cbind, mean_design_list)

  if (is.null(colnames(W_mean))) {
    colnames(W_mean) <- paste0("mean", seq_len(ncol(W_mean)))
  }

  colnames(W_mean) <- make.unique(colnames(W_mean))

  mean_groups <- list()
  mean_offset <- 0L

  for (expert_name in active_experts) {
    width_now <- ncol(mean_design_list[[expert_name]])
    mean_groups[[expert_name]] <- mean_offset + seq_len(width_now)
    mean_offset <- mean_offset + width_now
  }

  Z_mean <- cbind(beta0 = 1, W_mean)
  colnames(Z_mean) <- c(
    "beta0",
    paste0("gamma_", colnames(W_mean))
  )

  p_mean <- ncol(W_mean)
  p_coef <- ncol(Z_mean)
  prior_mean_coef <- c(m_beta, rep(0, p_mean))
  names(prior_mean_coef) <- colnames(Z_mean)

  theta_index <- list()
  theta_names <- character(0)
  theta_offset <- 0L

  for (expert_name in active_experts) {
    theta_index[[expert_name]] <- theta_offset + 1:2
    theta_names <- c(
      theta_names,
      paste0("log_ell_", expert_name),
      paste0("log_sigma_", expert_name)
    )
    theta_offset <- theta_offset + 2L
  }

  theta_dimension <- 2L * n_experts

  D2_cov <- if ("cov" %in% active_experts) {
    sqdist(X_cov)
  } else {
    NULL
  }

  network_metric_for_bounds <- if ("network" %in% active_experts) {
    network_metric_fun_factory(network_geometry)
  } else {
    NULL
  }

  hyperpriors <- prior_bundle$active_hyperpriors
  scale_tau <- prior_bundle$tau_scale[active_experts]
  tau2_default <- prior_bundle$tau2_default[active_experts]

  if (
    any(!is.finite(scale_tau)) ||
    any(scale_tau <= 0) ||
    any(!is.finite(tau2_default)) ||
    any(tau2_default <= 0)
  ) {
    stop("The active tau prior settings are invalid.")
  }

  prior_sd_blocks <- list()
  prior_sd_theta <- numeric(theta_dimension)

  for (expert_name in active_experts) {
    hp <- hyperpriors[[expert_name]]

    if (is.null(hp)) {
      stop("Missing hyperprior for expert ", expert_name, ".")
    }

    prior_sd_blocks[[expert_name]] <- pmax(
      c(hp$log_ell_sd, hp$log_sig_sd),
      1e-8
    )

    prior_sd_theta[theta_index[[expert_name]]] <-
      prior_sd_blocks[[expert_name]]
  }

  names(prior_sd_theta) <- theta_names

  ## ------------------------------------------------------------------------
  ## Expert kernels, bounds and priors
  ## ------------------------------------------------------------------------

  kernel_from_block <- function(expert_name, theta_block) {
    if (identical(expert_name, "cov")) {
      return(rbf_kernel_from_D2(
        D2_cov,
        sigma_f = exp(theta_block[2L]),
        l = exp(theta_block[1L])
      ))
    }

    network_kernel_from_geom(
      network_geometry,
      ell = exp(theta_block[1L]),
      sigma = exp(theta_block[2L])
    )
  }

  theta_ok_block <- function(expert_name, theta_block) {
    theta_block <- as.numeric(theta_block)

    if (length(theta_block) != 2L || any(!is.finite(theta_block))) {
      return(FALSE)
    }

    if (identical(expert_name, "cov")) {
      return(
        theta_block[1L] >= logl_cov_bounds[1L] &&
          theta_block[1L] <= logl_cov_bounds[2L] &&
          theta_block[2L] >= min_log_sigma &&
          theta_block[2L] <= max_log_sigma_cov
      )
    }

    if (
      theta_block[1L] < logl_network_bounds[1L] ||
      theta_block[1L] > logl_network_bounds[2L] ||
      theta_block[2L] < min_log_sigma ||
      theta_block[2L] > max_log_sigma_network
    ) {
      return(FALSE)
    }

    unit_metric <- network_metric_for_bounds(
      exp(theta_block[1L])
    )

    mean_marginal_sd <-
      exp(theta_block[2L]) * unit_metric$mean_sd

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

    pi_spike_now <- min(max(as.numeric(pi_spike), 0), 1)

    if (pi_spike_now <= 0) return(log_slab)

    spike_sd <- max(as.numeric(spike_log_sig_sd), 1e-8)

    log_spike <- dnorm(
      log_sigma,
      mean = spike_log_sig_mean,
      sd = spike_sd,
      log = TRUE
    )

    if (pi_spike_now >= 1) return(log_spike)

    aa <- log(pi_spike_now) + log_spike
    bb <- log(1 - pi_spike_now) + log_slab
    mm <- max(aa, bb)
    mm + log(exp(aa - mm) + exp(bb - mm))
  }

  logprior_block <- function(expert_name, theta_block) {
    if (!theta_ok_block(expert_name, theta_block)) return(-Inf)

    hp <- hyperpriors[[expert_name]]

    dnorm(
      theta_block[1L],
      mean = hp$log_ell_mean,
      sd = max(hp$log_ell_sd, 1e-8),
      log = TRUE
    ) +
      spike_slab_logsigma(theta_block[2L], hp)
  }

  logprior_all <- function(theta_vector) {
    theta_vector <- as.numeric(theta_vector)

    if (length(theta_vector) != theta_dimension) {
      return(-Inf)
    }

    value <- 0

    for (expert_name in active_experts) {
      value <- value + logprior_block(
        expert_name,
        theta_vector[theta_index[[expert_name]]]
      )
    }

    if (!is.finite(value)) return(-Inf)
    value
  }

  make_initial_theta_block <- function(expert_name) {
    hp <- hyperpriors[[expert_name]]
    theta <- c(hp$log_ell_mean, hp$log_sig_mean)

    if (identical(expert_name, "cov")) {
      theta[1L] <- min(
        max(theta[1L], logl_cov_bounds[1L]),
        logl_cov_bounds[2L]
      )
      theta[2L] <- min(
        max(theta[2L], min_log_sigma),
        max_log_sigma_cov
      )
    } else {
      theta[1L] <- min(
        max(theta[1L], logl_network_bounds[1L]),
        logl_network_bounds[2L]
      )
      theta[2L] <- min(
        max(theta[2L], min_log_sigma),
        max_log_sigma_network
      )

      while (
        !theta_ok_block(expert_name, theta) &&
        theta[2L] > min_log_sigma
      ) {
        theta[2L] <- max(
          theta[2L] - log(1.25),
          min_log_sigma
        )
      }
    }

    if (!theta_ok_block(expert_name, theta)) {
      stop("Initial theta is invalid for ", expert_name, ".")
    }

    theta
  }

  Csum_svd <- function(theta_vector) {
    K_sum <- matrix(0, n, n)

    for (expert_name in active_experts) {
      idx <- theta_index[[expert_name]]
      K_sum <- K_sum + kernel_from_block(
        expert_name,
        theta_vector[idx]
      )
    }

    K_sum <- symmetrize(K_sum) + diag(jitter, n)
    svd_psd_fallback(K_sum)
  }

  ## ------------------------------------------------------------------------
  ## Initial state
  ## ------------------------------------------------------------------------

  theta_vector <- numeric(theta_dimension)

  for (expert_name in active_experts) {
    theta_vector[theta_index[[expert_name]]] <-
      make_initial_theta_block(expert_name)
  }

  names(theta_vector) <- theta_names

  coef_vec <- prior_mean_coef
  tau2_group <- tau2_default
  names(tau2_group) <- active_experts

  aux_tau <- setNames(numeric(n_experts), active_experts)

  for (expert_name in active_experts) {
    aux_tau[expert_name] <- rInvGamma(
      1L,
      shape = 1 / 2,
      scale = 1 / (scale_tau[expert_name]^2)
    )
  }

  eta <- if (use_pu) {
    pmin(pmax(initial_eta, ETA_EPS), 1 - ETA_EPS)
  } else {
    0
  }

  mu_vec <- as.vector(Z_mean %*% coef_vec)
  f <- mu_vec

  svd_sum <- Csum_svd(theta_vector)
  U_sum <- svd_sum$U
  gam_sum <- pmax(svd_sum$d, 0)

  likelihood_log <- function(f_now, eta_now) {
    if (use_pu) {
      pu_loglik_vec(f_now, Y, eta_now)
    } else {
      bernoulli_loglik_vec(f_now, Y)
    }
  }

  likelihood_grad <- function(f_now, eta_now) {
    if (use_pu) {
      pu_gradloglik_vec(f_now, Y, eta_now)
    } else {
      bernoulli_gradloglik_vec(f_now, Y)
    }
  }

  log_likelihood <- likelihood_log(f, eta)
  gradient_f <- likelihood_grad(f, eta)

  kappa <- setNames(rep(as.numeric(kappa_init), n_experts), active_experts)
  proposal_cor <- setNames(
    lapply(seq_len(n_experts), function(ii) diag(2)),
    active_experts
  )

  joint_theta_scale <- 1 / sqrt(n_experts)
  delta_state <- as.numeric(delta_init)

  ## ------------------------------------------------------------------------
  ## Conditional coefficient and tau updates
  ## ------------------------------------------------------------------------

  prior_var_from_tau <- function(tau2_current) {
    output <- numeric(p_mean)

    for (expert_name in active_experts) {
      output[mean_groups[[expert_name]]] <- tau2_current[expert_name]
    }

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
      crossprod(Z_mean, KinvZ) +
        diag(1 / prior_var, p_coef)
    )

    information <- as.vector(crossprod(Z_mean, Kinvf)) +
      prior_mean_coef / prior_var

    R <- chol_safe(
      precision,
      jitter = 1e-10,
      max_tries = 8L
    )$R

    mean_coef <- as.vector(
      backsolve(
        R,
        forwardsolve(t(R), information)
      )
    )

    mean_coef + as.vector(backsolve(R, rnorm(p_coef)))
  }

  update_coef_tau <- function(
      f,
      U,
      gam,
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

      for (expert_name in active_experts) {
        idx <- mean_groups[[expert_name]]

        tau2_current[expert_name] <- rInvGamma(
          1L,
          shape = (nu_tau + length(idx)) / 2,
          scale =
            (nu_tau / aux_current[expert_name]) +
              sum(gamma_new[idx]^2) / 2
        )

        aux_current[expert_name] <- rInvGamma(
          1L,
          shape = (nu_tau + 1) / 2,
          scale =
            (1 / scale_tau[expert_name]^2) +
              (nu_tau / tau2_current[expert_name])
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

  update_remaining <- function(
      f,
      eta_current,
      U,
      gam,
      tau2_current,
      aux_current
  ) {
    eta_new <- eta_current
    eta_moved <- 0
    eta_n_eval <- 0L

    if (use_pu) {
      eta_update <- slice_eta_step(
        eta_current,
        f,
        Y,
        a_eta,
        b_eta,
        width = eta_slice_width,
        m = 100L
      )

      eta_new <- eta_update$eta
      eta_moved <- eta_update$moved
      eta_n_eval <- eta_update$n_eval
    }

    coef_update <- update_coef_tau(
      f,
      U,
      gam,
      tau2_current,
      aux_current
    )

    list(
      eta = eta_new,
      log_likelihood = likelihood_log(f, eta_new),
      gradient_f = likelihood_grad(f, eta_new),
      coef = coef_update$coef,
      tau2 = coef_update$tau2,
      aux_tau = coef_update$aux_tau,
      mu = coef_update$mu,
      eta_moved = eta_moved,
      eta_n_eval = eta_n_eval
    )
  }

  ## ------------------------------------------------------------------------
  ## f, blocked-theta and joint updates
  ## ------------------------------------------------------------------------

  agrad_f_step <- function(
      f,
      mu,
      log_likelihood,
      gradient_f,
      U,
      gam,
      eta,
      delta
  ) {
    nn <- length(f)
    residual <- f - mu

    z <- residual +
      (delta / 2) * gradient_f +
      sqrt(delta / 2) * rnorm(nn)

    sqrt_lambda_tilde <- sqrt(
      (pmax(gam, 0) * delta) /
        pmax(
          delta + 2 * pmax(gam, 0),
          .Machine$double.eps
        )
    )

    sqrt_lambda_tilde[!is.finite(sqrt_lambda_tilde)] <- 0

    temp1 <- crossprod(U, (2 / delta) * z)
    temp2 <- sqrt_lambda_tilde * temp1 + rnorm(nn)
    temp3 <- sqrt_lambda_tilde * temp2
    residual_proposed <- as.vector(U %*% temp3)
    f_proposed <- mu + residual_proposed

    log_likelihood_proposed <- likelihood_log(f_proposed, eta)
    gradient_proposed <- likelihood_grad(f_proposed, eta)

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
      expert_name,
      theta_current,
      f_fixed,
      mu_fixed,
      U_current,
      gam_current,
      kappa_scalar,
      correlation_block
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
      prior_sd = prior_sd_blocks[[expert_name]],
      cor_mat = correlation_block
    )

    proposal_increment <- rmvnorm0_chol(proposal_covariance)
    theta_proposed <- theta_current
    idx <- theta_index[[expert_name]]
    theta_proposed[idx] <- theta_proposed[idx] + proposal_increment

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
      eta,
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

    theta_proposed <- theta_current +
      rmvnorm0_chol(theta_covariance)

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
        pmax(
          delta + 2 * pmax(gam_proposed, 0),
          .Machine$double.eps
        )
    )

    sqrt_lambda_tilde[!is.finite(sqrt_lambda_tilde)] <- 0

    temp1 <- crossprod(U_proposed, (2 / delta) * z)
    temp2 <- sqrt_lambda_tilde * temp1 + rnorm(nn)
    temp3 <- sqrt_lambda_tilde * temp2
    residual_proposed <- as.vector(U_proposed %*% temp3)
    f_proposed <- mu + residual_proposed

    log_likelihood_proposed <- likelihood_log(f_proposed, eta)
    gradient_proposed <- likelihood_grad(f_proposed, eta)

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


  ## ------------------------------------------------------------------------
  ## Phase 1: exactly 10,000 correlated blocked warm-up iterations
  ## ------------------------------------------------------------------------

  start_time <- Sys.time()

  if (verbose) {
    cat("\n\n####################################################################\n")
    cat("# ", model_spec$method, "\n", sep = "")
    cat("####################################################################\n")
    cat("Active experts          :", paste(active_experts, collapse = " + "), "\n")
    cat("PU likelihood           :", use_pu, "\n")
    cat("Observations            :", n, "\n")
    cat("Mean-feature columns    :", p_mean, "\n")
    cat("Theta dimension         :", theta_dimension, "\n")
    cat("Blocked warm-up         :", n_blocked_burnin, "\n")
    cat("Joint warm-up           :", n_joint_burnin, "\n")
    cat("Posterior iterations    :", n_sample, "\n")
    cat("Thinning                :", thin, "\n")
    cat("Retained draws          :", n_sample %/% thin, "\n")
    cat("MCMC seed               :", model_spec$mcmc_seed, "\n")
  }

  acceptance_names <- c("f", active_experts)
  acceptance_blocked_sum <- setNames(
    numeric(length(acceptance_names)),
    acceptance_names
  )

  adaptation_window_length <- 0L
  adaptation_window_acceptance <- setNames(
    numeric(length(acceptance_names)),
    acceptance_names
  )
  adaptation_window_id <- 0L
  eta_moved_blocked_sum <- 0L

  theta_path_blocked <- matrix(
    NA_real_,
    nrow = n_blocked_burnin,
    ncol = theta_dimension,
    dimnames = list(NULL, theta_names)
  )

  cor_adapt_trace <- data.frame(
    iteration = integer(0),
    expert = character(0),
    rho = numeric(0),
    stringsAsFactors = FALSE
  )

  for (iteration in seq_len(n_blocked_burnin)) {
    f_update <- agrad_f_step(
      f,
      mu_vec,
      log_likelihood,
      gradient_f,
      U_sum,
      gam_sum,
      eta,
      delta_state
    )

    f <- f_update$f
    log_likelihood <- f_update$log_likelihood
    gradient_f <- f_update$gradient_f

    acceptance_now <- setNames(
      numeric(length(acceptance_names)),
      acceptance_names
    )
    acceptance_now["f"] <- as.integer(f_update$accepted)

    for (expert_name in active_experts) {
      theta_update <- theta_block_step(
        expert_name = expert_name,
        theta_current = theta_vector,
        f_fixed = f,
        mu_fixed = mu_vec,
        U_current = U_sum,
        gam_current = gam_sum,
        kappa_scalar = kappa[expert_name],
        correlation_block = proposal_cor[[expert_name]]
      )

      theta_vector <- theta_update$theta
      U_sum <- theta_update$U
      gam_sum <- theta_update$gam
      acceptance_now[expert_name] <-
        as.integer(theta_update$accepted)
    }

    remaining_update <- update_remaining(
      f,
      eta,
      U_sum,
      gam_sum,
      tau2_group,
      aux_tau
    )

    eta <- remaining_update$eta
    log_likelihood <- remaining_update$log_likelihood
    gradient_f <- remaining_update$gradient_f
    coef_vec <- remaining_update$coef
    tau2_group <- remaining_update$tau2
    aux_tau <- remaining_update$aux_tau
    mu_vec <- remaining_update$mu

    eta_moved_blocked_sum <- eta_moved_blocked_sum +
      remaining_update$eta_moved

    acceptance_blocked_sum <- acceptance_blocked_sum +
      acceptance_now

    theta_path_blocked[iteration, ] <- theta_vector

    if (
      adapt &&
      adapt_block_cor &&
      iteration >= cor_adapt_start &&
      iteration %% cor_adapt_interval == 0L
    ) {
      idx_cor <- seq.int(
        max(1L, iteration - cor_adapt_window + 1L),
        iteration
      )

      for (expert_name in active_experts) {
        idx_theta <- theta_index[[expert_name]]

        proposal_cor[[expert_name]] <- estimate_cor_2x2(
          theta_path_blocked[
            idx_cor,
            idx_theta,
            drop = FALSE
          ],
          fallback = proposal_cor[[expert_name]],
          shrinkage = cor_shrinkage,
          max_abs = cor_max_abs,
          min_n = min(200L, length(idx_cor))
        )

        cor_adapt_trace <- rbind(
          cor_adapt_trace,
          data.frame(
            iteration = iteration,
            expert = expert_name,
            rho = proposal_cor[[expert_name]][1L, 2L],
            stringsAsFactors = FALSE
          )
        )
      }
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

        for (expert_name in active_experts) {
          kappa[expert_name] <- adapt_log_step(
            kappa[expert_name],
            acceptance_rate[expert_name],
            target_theta,
            adaptation_window_id,
            adapt_rate_theta,
            t0_adapt,
            power_adapt,
            kappa_min,
            kappa_max
          )
        }

        adaptation_window_length <- 0L
        adaptation_window_acceptance[] <- 0
      }
    }

    if (
      verbose &&
      blocked_progress_every > 0L &&
      iteration %% blocked_progress_every == 0L
    ) {
      progress_message(
        paste0(model_spec$method, " blocked warm-up"),
        iteration,
        n_blocked_burnin
      )
    }
  }

  acceptance_blocked <-
    acceptance_blocked_sum / n_blocked_burnin

  eta_moved_blocked <- if (use_pu) {
    eta_moved_blocked_sum / n_blocked_burnin
  } else {
    NA_real_
  }

  proposal_cor_final <- lapply(
    proposal_cor,
    sanitize_cor_2x2,
    max_abs = cor_max_abs
  )

  ## ------------------------------------------------------------------------
  ## Construct the full joint proposal geometry from the blocked tail
  ## ------------------------------------------------------------------------

  fallback_correlation <- diag(theta_dimension)

  for (expert_name in active_experts) {
    idx <- theta_index[[expert_name]]
    fallback_correlation[idx, idx] <-
      proposal_cor_final[[expert_name]]
  }

  dimnames(fallback_correlation) <- list(theta_names, theta_names)

  idx_start <- max(
    1L,
    n_blocked_burnin - emp_cov_tail + 1L
  )

  theta_tail_blocked <- theta_path_blocked[
    idx_start:n_blocked_burnin,
    ,
    drop = FALSE
  ]

  joint_correlation_base <- estimate_correlation_matrix(
    theta_tail_blocked,
    fallback = fallback_correlation,
    shrinkage = joint_cor_shrinkage,
    max_abs = joint_cor_max_abs,
    min_n = min(200L, nrow(theta_tail_blocked))
  )

  block_scale <- numeric(theta_dimension)

  for (expert_name in active_experts) {
    block_scale[theta_index[[expert_name]]] <-
      sqrt(kappa[expert_name])
  }

  base_sd <- prior_sd_theta * block_scale
  D_theta <- diag(base_sd, theta_dimension)

  Sigma_theta_base <- D_theta %*%
    joint_correlation_base %*%
    D_theta

  Sigma_theta_base <- sanitize_covariance_matrix(
    Sigma_theta_base +
      cov_jitter * diag(theta_dimension)
  )

  dimnames(joint_correlation_base) <- list(theta_names, theta_names)
  dimnames(Sigma_theta_base) <- list(theta_names, theta_names)

  delta_joint <- delta_state

  if (verbose) {
    cat(sprintf(
      "[%s] %s blocked phase complete.\n",
      timestamp_now(),
      model_spec$method
    ))
    cat("  delta_joint          :", delta_joint, "\n")
    cat("  kappa                :\n")
    print(kappa)
    cat("  blocked acceptance   :\n")
    print(acceptance_blocked)
  }

  ## ------------------------------------------------------------------------
  ## Phase 2: exactly 20,000 joint warm-up/tuning iterations
  ## ------------------------------------------------------------------------

  acceptance_joint_tune_sum <- 0L
  acceptance_f_refresh_tune_sum <- 0L
  eta_moved_joint_tune_sum <- 0L
  joint_window_length <- 0L
  joint_window_acceptance <- 0L
  joint_window_id <- 0L

  for (iteration in seq_len(n_joint_burnin)) {
    Sigma_theta_effective <-
      (joint_theta_scale^2) * Sigma_theta_base

    joint_update <- joint_step(
      f,
      mu_vec,
      log_likelihood,
      gradient_f,
      theta_vector,
      eta,
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
        eta,
        delta_joint
      )

      f <- f_update$f
      log_likelihood <- f_update$log_likelihood
      gradient_f <- f_update$gradient_f

      acceptance_f_refresh_tune_sum <-
        acceptance_f_refresh_tune_sum +
        as.integer(f_update$accepted)
    }

    remaining_update <- update_remaining(
      f,
      eta,
      U_sum,
      gam_sum,
      tau2_group,
      aux_tau
    )

    eta <- remaining_update$eta
    log_likelihood <- remaining_update$log_likelihood
    gradient_f <- remaining_update$gradient_f
    coef_vec <- remaining_update$coef
    tau2_group <- remaining_update$tau2
    aux_tau <- remaining_update$aux_tau
    mu_vec <- remaining_update$mu

    eta_moved_joint_tune_sum <-
      eta_moved_joint_tune_sum +
      remaining_update$eta_moved

    if (adapt && iteration >= adapt_start) {
      joint_window_length <- joint_window_length + 1L
      joint_window_acceptance <-
        joint_window_acceptance + accepted_joint

      if (joint_window_length == adapt_interval) {
        joint_window_id <- joint_window_id + 1L
        joint_rate <-
          joint_window_acceptance / joint_window_length

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

    if (
      verbose &&
      joint_progress_every > 0L &&
      iteration %% joint_progress_every == 0L
    ) {
      progress_message(
        paste0(model_spec$method, " joint warm-up"),
        iteration,
        n_joint_burnin
      )
    }
  }

  acceptance_joint_tune <-
    acceptance_joint_tune_sum / n_joint_burnin

  acceptance_f_refresh_tune <-
    acceptance_f_refresh_tune_sum / n_joint_burnin

  eta_moved_joint_tune <- if (use_pu) {
    eta_moved_joint_tune_sum / n_joint_burnin
  } else {
    NA_real_
  }

  delta_final <- delta_joint
  joint_theta_scale_final <- joint_theta_scale

  Sigma_theta_final <- sanitize_covariance_matrix(
    (joint_theta_scale_final^2) * Sigma_theta_base
  )

  dimnames(Sigma_theta_final) <- list(theta_names, theta_names)

  block_covariance_final <- list()

  for (expert_name in active_experts) {
    block_covariance_final[[expert_name]] <-
      sanitize_covariance_matrix(
        proposal_cov_scaled_2x2(
          kappa[expert_name],
          prior_sd_blocks[[expert_name]],
          proposal_cor_final[[expert_name]],
          jitter = cov_jitter
        )
      )
  }

  ## ------------------------------------------------------------------------
  ## Phase 3: 150,000 frozen posterior iterations, thinning 30
  ## ------------------------------------------------------------------------

  n_saved_expected <- n_sample %/% thin
  n_saved <- 0L

  latent_probability_draws <- matrix(
    NA_real_,
    nrow = n_saved_expected,
    ncol = n,
    dimnames = list(
      NULL,
      paste0("node_", seq_len(n))
    )
  )

  observed_probability_draws <- matrix(
    NA_real_,
    nrow = n_saved_expected,
    ncol = n,
    dimnames = dimnames(latent_probability_draws)
  )

  conditional_T_probability_draws <- matrix(
    NA_real_,
    nrow = n_saved_expected,
    ncol = n,
    dimnames = dimnames(latent_probability_draws)
  )

  eta_draws <- if (use_pu) {
    numeric(n_saved_expected)
  } else {
    rep(NA_real_, n_saved_expected)
  }

  beta0_draws <- numeric(n_saved_expected)

  tau_draws <- matrix(
    NA_real_,
    nrow = n_saved_expected,
    ncol = n_experts,
    dimnames = list(
      NULL,
      paste0("tau_", active_experts)
    )
  )

  theta_draw_names <- unlist(
    lapply(active_experts, function(expert_name) {
      c(
        paste0("ell_", expert_name),
        paste0("sigma_", expert_name)
      )
    }),
    use.names = FALSE
  )

  theta_draws <- matrix(
    NA_real_,
    nrow = n_saved_expected,
    ncol = theta_dimension,
    dimnames = list(NULL, theta_draw_names)
  )

  acceptance_joint_sampling_sum <- 0L
  acceptance_f_refresh_sampling_sum <- 0L
  eta_moved_sampling_sum <- 0L

  zero_idx <- which(Y == 0L)
  observed_positive_idx <- which(Y == 1L)

  for (iteration in seq_len(n_sample)) {
    joint_update <- joint_step(
      f,
      mu_vec,
      log_likelihood,
      gradient_f,
      theta_vector,
      eta,
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
        eta,
        delta_final
      )

      f <- f_update$f
      log_likelihood <- f_update$log_likelihood
      gradient_f <- f_update$gradient_f

      acceptance_f_refresh_sampling_sum <-
        acceptance_f_refresh_sampling_sum +
        as.integer(f_update$accepted)
    }

    remaining_update <- update_remaining(
      f,
      eta,
      U_sum,
      gam_sum,
      tau2_group,
      aux_tau
    )

    eta <- remaining_update$eta
    log_likelihood <- remaining_update$log_likelihood
    gradient_f <- remaining_update$gradient_f
    coef_vec <- remaining_update$coef
    tau2_group <- remaining_update$tau2
    aux_tau <- remaining_update$aux_tau
    mu_vec <- remaining_update$mu

    eta_moved_sampling_sum <-
      eta_moved_sampling_sum +
      remaining_update$eta_moved

    if (iteration %% thin == 0L) {
      n_saved <- n_saved + 1L
      p_latent <- plogis(f)

      latent_probability_draws[n_saved, ] <- p_latent

      if (use_pu) {
        observed_probability_draws[n_saved, ] <-
          (1 - eta) * p_latent

        conditional_T_probability_draws[
          n_saved,
          zero_idx
        ] <- hidden_probability_given_zero(
          p_latent[zero_idx],
          eta
        )

        conditional_T_probability_draws[
          n_saved,
          observed_positive_idx
        ] <- 1

        eta_draws[n_saved] <- eta
      } else {
        observed_probability_draws[n_saved, ] <- p_latent
        conditional_T_probability_draws[n_saved, ] <- p_latent
      }

      beta0_draws[n_saved] <- coef_vec[1L]
      tau_draws[n_saved, ] <- sqrt(tau2_group[active_experts])
      theta_draws[n_saved, ] <- exp(theta_vector)
    }

    if (
      verbose &&
      sample_progress_every > 0L &&
      iteration %% sample_progress_every == 0L
    ) {
      progress_message(
        paste0(model_spec$method, " posterior sampling"),
        iteration,
        n_sample
      )
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

  if (
    any(!is.finite(latent_probability_draws)) ||
    any(!is.finite(observed_probability_draws)) ||
    any(!is.finite(conditional_T_probability_draws))
  ) {
    stop("At least one posterior probability matrix was not filled.")
  }

  end_time <- Sys.time()

  runtime_sec <- as.numeric(difftime(
    end_time,
    start_time,
    units = "secs"
  ))

  acceptance_joint_sampling <-
    acceptance_joint_sampling_sum / n_sample

  acceptance_f_refresh_sampling <-
    acceptance_f_refresh_sampling_sum / n_sample

  eta_moved_sampling <- if (use_pu) {
    eta_moved_sampling_sum / n_sample
  } else {
    NA_real_
  }

  if (verbose) {
    cat(sprintf(
      "[%s] %s finished.\n",
      timestamp_now(),
      model_spec$method
    ))
    cat("Runtime                    :", format_duration(runtime_sec), "\n")
    cat("Joint tune acceptance      :", acceptance_joint_tune, "\n")
    cat("Joint sample acceptance    :", acceptance_joint_sampling, "\n")
    cat("f-refresh sample accept.   :", acceptance_f_refresh_sampling, "\n")
    if (use_pu) {
      cat("eta move rate              :", eta_moved_sampling, "\n")
      cat("Final eta                  :", eta, "\n")
    }
    cat("Final theta (ell/sigma)    :\n")
    print(setNames(exp(theta_vector), theta_draw_names))
    cat("Final tau                 :\n")
    print(sqrt(tau2_group))
    cat("Final delta               :", delta_final, "\n")
    cat("Final joint theta scale   :", joint_theta_scale_final, "\n")
  }

  list(
    model_spec = model_spec,
    state = list(
      theta_log = setNames(as.numeric(theta_vector), theta_names),
      beta = setNames(as.numeric(coef_vec), colnames(Z_mean)),
      tau2 = setNames(as.numeric(tau2_group), active_experts),
      aux_tau = setNames(as.numeric(aux_tau), active_experts),
      eta = if (use_pu) as.numeric(eta) else NA_real_,
      f = as.numeric(f)
    ),
    draws = list(
      latent_probability = latent_probability_draws,
      observed_probability = observed_probability_draws,
      conditional_T_probability = conditional_T_probability_draws,
      eta = eta_draws,
      beta0 = beta0_draws,
      tau = tau_draws,
      theta = theta_draws
    ),
    proposal = list(
      delta_final = as.numeric(delta_final),
      kappa_final = kappa,
      proposal_cor_final = proposal_cor_final,
      block_covariance_final = block_covariance_final,
      joint_correlation_base = joint_correlation_base,
      Sigma_theta_base = Sigma_theta_base,
      joint_theta_scale_final = as.numeric(joint_theta_scale_final),
      Sigma_theta_final = Sigma_theta_final,
      cor_adapt_trace = cor_adapt_trace
    ),
    acceptance = list(
      blocked = acceptance_blocked,
      joint_tune = acceptance_joint_tune,
      f_refresh_tune = acceptance_f_refresh_tune,
      joint_sampling = acceptance_joint_sampling,
      f_refresh_sampling = acceptance_f_refresh_sampling,
      eta_moved_blocked = eta_moved_blocked,
      eta_moved_joint_tune = eta_moved_joint_tune,
      eta_moved_sampling = eta_moved_sampling
    ),
    phases = c(
      blocked_warmup = n_blocked_burnin,
      joint_warmup = n_joint_burnin,
      posterior_iterations = n_sample,
      thinning = thin,
      retained = n_saved
    ),
    time_start = start_time,
    time_end = end_time,
    runtime_sec = runtime_sec
  )
}


## ============================================================================
## 11. HPD SUMMARIES FROM CODA
## ============================================================================

hpd_interval_coda <- function(x, probability = HPD_PROBABILITY) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    stop("At least two finite draws are required for an HPD interval.")
  }

  if (!is.finite(probability) || probability <= 0 || probability >= 1) {
    stop("The HPD probability must lie strictly between zero and one.")
  }

  ## coda::HPDinterval may receive a constant chain (e.g. Pr(T=1|Y=1)=1 in
  ## the PU model). In that special case the exact HPD interval is degenerate.
  if (max(x) - min(x) <= 10 * .Machine$double.eps) {
    value <- x[1L]
    return(c(lower = value, upper = value))
  }

  interval <- coda::HPDinterval(
    coda::as.mcmc(x),
    prob = probability
  )

  c(
    lower = as.numeric(interval[1L, "lower"]),
    upper = as.numeric(interval[1L, "upper"])
  )
}

hpd_matrix_coda <- function(
    draw_matrix,
    probability = HPD_PROBABILITY
) {
  draw_matrix <- as.matrix(draw_matrix)

  if (
    nrow(draw_matrix) < 2L ||
    ncol(draw_matrix) < 1L ||
    any(!is.finite(draw_matrix))
  ) {
    stop("The draw matrix supplied to coda HPD calculation is invalid.")
  }

  intervals <- t(vapply(
    seq_len(ncol(draw_matrix)),
    function(column_index) {
      hpd_interval_coda(
        draw_matrix[, column_index],
        probability = probability
      )
    },
    numeric(2)
  ))

  colnames(intervals) <- c("lower", "upper")
  intervals
}

roc_auc <- function(scores, labels) {
  scores <- as.numeric(scores)
  labels <- as.integer(labels)
  ok <- is.finite(scores) & labels %in% c(0L, 1L)
  scores <- scores[ok]
  labels <- labels[ok]

  n_positive <- sum(labels == 1L)
  n_negative <- sum(labels == 0L)

  if (n_positive == 0L || n_negative == 0L) {
    return(NA_real_)
  }

  ranks_ascending <- rank(scores, ties.method = "average")
  sum_ranks_positive <- sum(ranks_ascending[labels == 1L])

  (
    sum_ranks_positive -
      n_positive * (n_positive + 1) / 2
  ) / (
    n_positive * n_negative
  )
}

overall_logloss <- function(
    probabilities,
    truth,
    eps = 1e-8
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
  mean(
    (as.numeric(probabilities) - as.numeric(truth))^2
  )
}

summarize_completed_model <- function(
    fit,
    model_spec,
    prior_bundle
) {
  latent_draws <- as.matrix(fit$draws$latent_probability)
  observed_draws <- as.matrix(fit$draws$observed_probability)
  conditional_T_draws <- as.matrix(
    fit$draws$conditional_T_probability
  )

  expected_dimension <- c(N_RETAINED, N_NODES)

  if (
    !all(dim(latent_draws) == expected_dimension) ||
    !all(dim(observed_draws) == expected_dimension) ||
    !all(dim(conditional_T_draws) == expected_dimension)
  ) {
    stop(
      model_spec$method,
      " did not return the expected 5,000 x 100 probability matrices."
    )
  }

  latent_mean <- colMeans(latent_draws)
  observed_mean <- colMeans(observed_draws)
  conditional_T_mean <- colMeans(conditional_T_draws)

  latent_hpd <- hpd_matrix_coda(
    latent_draws,
    probability = HPD_PROBABILITY
  )

  observed_hpd <- hpd_matrix_coda(
    observed_draws,
    probability = HPD_PROBABILITY
  )

  conditional_T_hpd <- hpd_matrix_coda(
    conditional_T_draws,
    probability = HPD_PROBABILITY
  )

  node_summary <- data.frame(
    model_id = model_spec$model_id,
    method = model_spec$method,
    active_experts = paste(
      model_spec$active_experts,
      collapse = " + "
    ),
    use_PU_likelihood = isTRUE(model_spec$use_pu),
    node = seq_len(N_NODES),
    T = T_TRUE,
    Y = Y,
    node_group = NODE_GROUP,

    ## Primary figure-ready posterior summaries: p_i = plogis(f_i).
    posterior_predictive_mean = latent_mean,
    HPD80_lower = latent_hpd[, "lower"],
    HPD80_upper = latent_hpd[, "upper"],
    HPD80_width = latent_hpd[, "upper"] - latent_hpd[, "lower"],

    latent_probability_mean = latent_mean,
    latent_HPD80_lower = latent_hpd[, "lower"],
    latent_HPD80_upper = latent_hpd[, "upper"],
    latent_HPD80_width = latent_hpd[, "upper"] - latent_hpd[, "lower"],

    ## Probability of the observed label under each fitted model.
    observed_probability_mean = observed_mean,
    observed_HPD80_lower = observed_hpd[, "lower"],
    observed_HPD80_upper = observed_hpd[, "upper"],
    observed_HPD80_width =
      observed_hpd[, "upper"] - observed_hpd[, "lower"],

    ## Draw-by-draw Pr(T_i=1 | Y_i, parameters). This differs from p_i only
    ## for GP-PU-PoE.
    conditional_T_probability_mean = conditional_T_mean,
    conditional_T_HPD80_lower = conditional_T_hpd[, "lower"],
    conditional_T_HPD80_upper = conditional_T_hpd[, "upper"],
    conditional_T_HPD80_width =
      conditional_T_hpd[, "upper"] - conditional_T_hpd[, "lower"],

    stringsAsFactors = FALSE
  )

  unlabeled_idx <- which(Y == 0L)
  hidden_truth_unlabeled <- as.integer(T_TRUE[unlabeled_idx] == 1L)

  eta_summary <- if (isTRUE(model_spec$use_pu)) {
    eta_hpd <- hpd_interval_coda(
      fit$draws$eta,
      probability = HPD_PROBABILITY
    )

    data.frame(
      model_id = model_spec$model_id,
      method = model_spec$method,
      eta_true = length(HIDDEN_POS_IDX) / length(TRUE_POS_IDX),
      eta_posterior_mean = mean(fit$draws$eta),
      eta_posterior_sd = stats::sd(fit$draws$eta),
      eta_HPD80_lower = eta_hpd["lower"],
      eta_HPD80_upper = eta_hpd["upper"],
      eta_HPD80_width = eta_hpd["upper"] - eta_hpd["lower"],
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      model_id = model_spec$model_id,
      method = model_spec$method,
      eta_true = length(HIDDEN_POS_IDX) / length(TRUE_POS_IDX),
      eta_posterior_mean = NA_real_,
      eta_posterior_sd = NA_real_,
      eta_HPD80_lower = NA_real_,
      eta_HPD80_upper = NA_real_,
      eta_HPD80_width = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  blocked_acceptance <- fit$acceptance$blocked

  diagnostics <- data.frame(
    model_id = model_spec$model_id,
    method = model_spec$method,
    active_experts = paste(
      model_spec$active_experts,
      collapse = " + "
    ),
    use_PU_likelihood = isTRUE(model_spec$use_pu),
    n_nodes = N_NODES,
    n_covariates = N_COVARIATES,
    n_network_features = ncol(Z_network),
    blocked_burnin = fit$phases["blocked_warmup"],
    joint_burnin = fit$phases["joint_warmup"],
    posterior_iterations = fit$phases["posterior_iterations"],
    thinning = fit$phases["thinning"],
    retained_draws = fit$phases["retained"],
    acc_blocked_f = unname(blocked_acceptance["f"]),
    acc_blocked_cov = if ("cov" %in% names(blocked_acceptance)) {
      unname(blocked_acceptance["cov"])
    } else {
      NA_real_
    },
    acc_blocked_network = if ("network" %in% names(blocked_acceptance)) {
      unname(blocked_acceptance["network"])
    } else {
      NA_real_
    },
    acc_joint_tune = fit$acceptance$joint_tune,
    acc_f_refresh_tune = fit$acceptance$f_refresh_tune,
    acc_joint_sampling = fit$acceptance$joint_sampling,
    acc_f_refresh_sampling = fit$acceptance$f_refresh_sampling,
    eta_move_rate_sampling = fit$acceptance$eta_moved_sampling,
    delta_final = fit$proposal$delta_final,
    joint_theta_scale_final = fit$proposal$joint_theta_scale_final,
    runtime_sec = fit$runtime_sec,
    overall_AUC_latent_mean = roc_auc(latent_mean, T_TRUE),
    hidden_AUC_latent_mean = roc_auc(
      latent_mean[unlabeled_idx],
      hidden_truth_unlabeled
    ),
    overall_LogLoss_latent_mean = overall_logloss(
      latent_mean,
      T_TRUE
    ),
    overall_Brier_latent_mean = overall_brier(
      latent_mean,
      T_TRUE
    ),
    mean_HPD80_width = mean(node_summary$HPD80_width),
    mean_HPD80_width_true_zero = mean(
      node_summary$HPD80_width[node_summary$node_group == "true zero"]
    ),
    mean_HPD80_width_hidden_positive = mean(
      node_summary$HPD80_width[
        node_summary$node_group == "hidden positive"
      ]
    ),
    mean_HPD80_width_observed_positive = mean(
      node_summary$HPD80_width[
        node_summary$node_group == "observed positive"
      ]
    ),
    stringsAsFactors = FALSE
  )

  list(
    node_summary = node_summary,
    eta_summary = eta_summary,
    diagnostics = diagnostics,
    prior_calibration = prior_bundle$summary
  )
}


## ============================================================================
## 12. RUN OR LOAD THE 11 DISTINCT PU-POE FITS
## ============================================================================

make_sensitivity_model_spec <- function(fit_spec) {
  list(
    model_id = fit_spec$fit_id,
    method = "GP-PU-PoE",
    active_experts = BASE_PUPOE_MODEL_SPEC$active_experts,
    use_pu = TRUE,
    mcmc_seed = as.integer(fit_spec$mcmc_seed)
  )
}

make_sensitivity_signature <- function(
    fit_spec,
    model_spec,
    prior_bundle
) {
  list(
    script_version = SCRIPT_VERSION,
    toy_data_hash = TOY_DATA_HASH,
    toy_seed = TOY_SEED,
    modularity_seed = MODULARITY_SEED,
    fit_spec = fit_spec,
    model_spec = model_spec,
    MCMC = list(
      blocked_burnin = BLOCKED_BURNIN,
      joint_burnin = JOINT_BURNIN,
      posterior_iterations = POSTERIOR_ITERATIONS,
      thin = THIN,
      retained = N_RETAINED
    ),
    data = list(
      n_nodes = N_NODES,
      n_covariates = N_COVARIATES,
      true_zero_idx = TRUE_ZERO_IDX,
      hidden_positive_idx = HIDDEN_POS_IDX,
      observed_positive_idx = OBSERVED_POS_IDX,
      p_00 = P_00,
      p_01 = P_01,
      p_11 = P_11,
      modularity_q = modularity$q,
      modularity_selected_eigenvalues =
        modularity$selected_eigenvalues
    ),
    priors = list(
      beta0 = prior_bundle$beta0_prior,
      eta = prior_bundle$eta_prior,
      active_hyperpriors = prior_bundle$active_hyperpriors,
      tau_scale = prior_bundle$tau_scale,
      effective_sd_center = prior_bundle$effective_sd_center,
      max_sigma_cov = fit_spec$max_sigma_cov,
      max_mean_sd_network = fit_spec$max_mean_sd_network
    ),
    requested_interval = list(
      type = "HPD",
      package = "coda",
      function_name = "HPDinterval",
      probability = HPD_PROBABILITY
    )
  )
}

add_fit_metadata <- function(df, fit_spec) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  data.frame(
    fit_id = fit_spec$fit_id,
    varying_block = fit_spec$varying_block,
    prior_role = fit_spec$prior_role,
    df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

run_or_load_sensitivity_fit <- function(fit_spec) {
  prior_bundle <- sensitivity_prior_bundles[[fit_spec$fit_id]]
  model_spec <- make_sensitivity_model_spec(fit_spec)
  fit_signature <- make_sensitivity_signature(
    fit_spec,
    model_spec,
    prior_bundle
  )

  checkpoint_file <- file.path(
    FIT_CHECKPOINT_DIR,
    paste0(fit_spec$fit_id, "_complete.rds")
  )

  if (
    isTRUE(RESUME_COMPLETED_FITS) &&
    !isTRUE(FORCE_FRESH_RUN) &&
    file.exists(checkpoint_file)
  ) {
    candidate <- tryCatch(
      readRDS(checkpoint_file),
      error = function(e) NULL
    )

    if (
      is.list(candidate) &&
      isTRUE(candidate$complete) &&
      identical(candidate$fit_signature, fit_signature)
    ) {
      cat(sprintf(
        "[%s] Reusing completed sensitivity fit %s.\n",
        timestamp_now(),
        fit_spec$fit_id
      ))
      return(candidate)
    }
  }

  cat("\n\n####################################################################\n")
  cat("# PRIOR SENSITIVITY FIT: ", fit_spec$fit_id, "\n", sep = "")
  cat("####################################################################\n")
  cat("Varying block          :", fit_spec$varying_block, "\n")
  cat("Prior role             :", fit_spec$prior_role, "\n")
  cat("Common MCMC seed       :", fit_spec$mcmc_seed, "\n")
  cat("Total GP SD centre     :", fit_spec$total_sd_center, "\n")
  cat("Correlation calibration:",
      paste(fit_spec$corr_pct_range, collapse = ", "), "\n")
  cat("Mean RMS multiplier    :", fit_spec$mean_rms_95_mult, "\n")
  cat("beta0 prevalence range :",
      paste(fit_spec$beta0_prev_range, collapse = ", "), "\n")
  cat("eta prior              : Beta(",
      prior_bundle$eta_prior$a_eta, ", ",
      prior_bundle$eta_prior$b_eta, ")\n", sep = "")

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

  summaries <- summarize_completed_model(
    fit = fit,
    model_spec = model_spec,
    prior_bundle = prior_bundle
  )

  result <- list(
    complete = TRUE,
    fit_signature = fit_signature,
    fit_spec = fit_spec,
    model_spec = model_spec,
    prior_bundle = prior_bundle,
    fit = fit,
    node_summary = add_fit_metadata(
      summaries$node_summary,
      fit_spec
    ),
    eta_summary = add_fit_metadata(
      summaries$eta_summary,
      fit_spec
    ),
    diagnostics = add_fit_metadata(
      summaries$diagnostics,
      fit_spec
    ),
    prior_calibration = add_fit_metadata(
      summaries$prior_calibration,
      fit_spec
    ),
    completed_at = timestamp_now()
  )

  saveRDS(result, checkpoint_file)
  result
}

sensitivity_results <- vector(
  "list",
  length(SENSITIVITY_RUN_ORDER)
)
names(sensitivity_results) <- SENSITIVITY_RUN_ORDER

all_fits_start <- Sys.time()

for (fit_id_now in SENSITIVITY_RUN_ORDER) {
  sensitivity_results[[fit_id_now]] <-
    run_or_load_sensitivity_fit(
      SENSITIVITY_SPECS[[fit_id_now]]
    )
}

all_fits_end <- Sys.time()

if (any(!vapply(
  sensitivity_results,
  function(x) is.list(x) && isTRUE(x$complete),
  logical(1)
))) {
  stop("Not all 11 prior-sensitivity fits completed successfully.")
}

## ============================================================================
## 13. COMBINE POSTERIOR SUMMARIES AND ROBUSTNESS DIAGNOSTICS
## ============================================================================

node_summaries_long <- do.call(
  rbind,
  lapply(sensitivity_results, `[[`, "node_summary")
)
row.names(node_summaries_long) <- NULL

eta_summary <- do.call(
  rbind,
  lapply(sensitivity_results, `[[`, "eta_summary")
)
row.names(eta_summary) <- NULL

## ============================================================================
## COMPACT ETA-PRIOR SENSITIVITY SUMMARY
## ============================================================================

ETA_PRIOR_SENSITIVITY_IDS <- c(
  "ETA_CONSERVATIVE",
  "DEFAULT",
  "ETA_DIFFUSE"
)

eta_prior_sensitivity_summary <- do.call(
  rbind,
  lapply(
    ETA_PRIOR_SENSITIVITY_IDS,
    function(fit_id_now) {
      
      eta_row <- eta_summary[
        eta_summary$fit_id == fit_id_now,
        ,
        drop = FALSE
      ]
      
      if (nrow(eta_row) != 1L) {
        stop(
          "Expected exactly one eta-summary row for ",
          fit_id_now,
          "."
        )
      }
      
      eta_prior_now <-
        sensitivity_results[[fit_id_now]]$prior_bundle$eta_prior
      
      data.frame(
        fit_id = fit_id_now,
        prior_role =
          sensitivity_results[[fit_id_now]]$fit_spec$prior_role,
        eta_prior_alpha =
          eta_prior_now$a_eta,
        eta_prior_beta =
          eta_prior_now$b_eta,
        eta_prior_mean =
          eta_prior_now$a_eta /
          (
            eta_prior_now$a_eta +
              eta_prior_now$b_eta
          ),
        eta_true =
          eta_row$eta_true,
        eta_posterior_mean =
          eta_row$eta_posterior_mean,
        eta_posterior_sd =
          eta_row$eta_posterior_sd,
        eta_HPD80_lower =
          eta_row$eta_HPD80_lower,
        eta_HPD80_upper =
          eta_row$eta_HPD80_upper,
        eta_HPD80_width =
          eta_row$eta_HPD80_width,
        stringsAsFactors = FALSE
      )
    }
  )
)

row.names(eta_prior_sensitivity_summary) <- NULL

write.csv(
  eta_prior_sensitivity_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_ETA_PRIOR_SENSITIVITY_summary_HPD80.csv"
  ),
  row.names = FALSE
)

cat("\n================ ETA-PRIOR SENSITIVITY ================\n")
print(
  eta_prior_sensitivity_summary,
  row.names = FALSE
)

for (
  row_index in seq_len(
    nrow(eta_prior_sensitivity_summary)
  )
) {
  row_now <- eta_prior_sensitivity_summary[
    row_index,
    ,
    drop = FALSE
  ]
  
  cat(sprintf(
    paste0(
      "%s: eta ~ Beta(%.2f, %.2f), prior mean = %.3f; ",
      "posterior mean = %.3f, 80%% HPD interval = [%.3f, %.3f], ",
      "true eta = %.3f.\n"
    ),
    row_now$prior_role,
    row_now$eta_prior_alpha,
    row_now$eta_prior_beta,
    row_now$eta_prior_mean,
    row_now$eta_posterior_mean,
    row_now$eta_HPD80_lower,
    row_now$eta_HPD80_upper,
    row_now$eta_true
  ))
}

fit_diagnostics <- do.call(
  rbind,
  lapply(sensitivity_results, `[[`, "diagnostics")
)
row.names(fit_diagnostics) <- NULL

prior_calibration_completed <- do.call(
  rbind,
  lapply(sensitivity_results, `[[`, "prior_calibration")
)
row.names(prior_calibration_completed) <- NULL

write.csv(
  node_summaries_long,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_node_summaries_HPD80.csv"
  ),
  row.names = FALSE
)

write.csv(
  eta_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_eta_HPD80.csv"
  ),
  row.names = FALSE
)

write.csv(
  fit_diagnostics,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_MCMC_diagnostics.csv"
  ),
  row.names = FALSE
)

write.csv(
  prior_calibration_completed,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_active_priors.csv"
  ),
  row.names = FALSE
)

baseline_nodes <- node_summaries_long[
  node_summaries_long$fit_id == "DEFAULT",
  c("node", "latent_probability_mean"),
  drop = FALSE
]
baseline_nodes <- baseline_nodes[order(baseline_nodes$node), , drop = FALSE]

if (nrow(baseline_nodes) != N_NODES) {
  stop("The default fit did not return exactly 100 node summaries.")
}

baseline_probability <- baseline_nodes$latent_probability_mean
names(baseline_probability) <- baseline_nodes$node

ordered_top_nodes <- function(probability, k) {
  node_index <- seq_along(probability)
  node_index[
    order(-probability, node_index)[seq_len(min(k, length(probability)))]
  ]
}

fit_summary_rows <- list()
topk_rows <- list()
fit_summary_index <- 1L
topk_index <- 1L

for (fit_id_now in SENSITIVITY_RUN_ORDER) {
  result_now <- sensitivity_results[[fit_id_now]]
  spec_now <- result_now$fit_spec
  prior_now <- result_now$prior_bundle

  node_now <- result_now$node_summary[
    order(result_now$node_summary$node),
    ,
    drop = FALSE
  ]
  probability_now <- node_now$latent_probability_mean

  spearman_now <- suppressWarnings(stats::cor(
    probability_now,
    baseline_probability,
    method = "spearman"
  ))

  eta_now <- result_now$eta_summary[1L, , drop = FALSE]
  diagnostics_now <- result_now$diagnostics[1L, , drop = FALSE]

  fit_summary_rows[[fit_summary_index]] <- data.frame(
    fit_id = fit_id_now,
    varying_block = spec_now$varying_block,
    prior_role = spec_now$prior_role,
    total_sd_center = spec_now$total_sd_center,
    corr_pct_low = spec_now$corr_pct_range[1L],
    corr_pct_high = spec_now$corr_pct_range[2L],
    mean_rms_95_mult = spec_now$mean_rms_95_mult,
    beta0_prev_low = spec_now$beta0_prev_range[1L],
    beta0_prev_high = spec_now$beta0_prev_range[2L],
    beta0_prior_mean = prior_now$beta0_prior$mean,
    beta0_prior_sd = prior_now$beta0_prior$sd,
    eta_prior_alpha = prior_now$eta_prior$a_eta,
    eta_prior_beta = prior_now$eta_prior$b_eta,
    eta_posterior_mean = eta_now$eta_posterior_mean,
    eta_posterior_sd = eta_now$eta_posterior_sd,
    eta_HPD80_lower = eta_now$eta_HPD80_lower,
    eta_HPD80_upper = eta_now$eta_HPD80_upper,
    spearman_with_default = spearman_now,
    mean_absolute_probability_difference = mean(
      abs(probability_now - baseline_probability)
    ),
    max_absolute_probability_difference = max(
      abs(probability_now - baseline_probability)
    ),
    mean_HPD80_width = mean(node_now$latent_HPD80_width),
    overall_AUC = diagnostics_now$overall_AUC_latent_mean,
    hidden_AUC = diagnostics_now$hidden_AUC_latent_mean,
    runtime_sec = result_now$fit$runtime_sec,
    stringsAsFactors = FALSE
  )
  fit_summary_index <- fit_summary_index + 1L

  for (k_now in TOP_K_VALUES) {
    current_top <- ordered_top_nodes(probability_now, k_now)
    baseline_top <- ordered_top_nodes(baseline_probability, k_now)
    overlap_count <- length(intersect(current_top, baseline_top))

    topk_rows[[topk_index]] <- data.frame(
      fit_id = fit_id_now,
      varying_block = spec_now$varying_block,
      prior_role = spec_now$prior_role,
      k = k_now,
      overlap_count = overlap_count,
      overlap_fraction = overlap_count / k_now,
      stringsAsFactors = FALSE
    )
    topk_index <- topk_index + 1L
  }
}

sensitivity_fit_summary <- do.call(rbind, fit_summary_rows)
row.names(sensitivity_fit_summary) <- NULL

topk_overlap_summary <- do.call(rbind, topk_rows)
row.names(topk_overlap_summary) <- NULL

summarize_parameter_draws <- function(result_now) {
  fit_id_now <- result_now$fit_spec$fit_id
  spec_now <- result_now$fit_spec
  draws_now <- result_now$fit$draws

  parameter_vectors <- list(
    beta0 = as.numeric(draws_now$beta0),
    eta = as.numeric(draws_now$eta),
    tau_cov = as.numeric(draws_now$tau[, "tau_cov"]),
    tau_network = as.numeric(draws_now$tau[, "tau_network"]),
    ell_cov = as.numeric(draws_now$theta[, "ell_cov"]),
    sigma_cov = as.numeric(draws_now$theta[, "sigma_cov"]),
    ell_network = as.numeric(draws_now$theta[, "ell_network"]),
    sigma_network = as.numeric(draws_now$theta[, "sigma_network"])
  )

  rows <- lapply(
    names(parameter_vectors),
    function(parameter_name) {
      values <- parameter_vectors[[parameter_name]]
      interval <- hpd_interval_coda(
        values,
        probability = HPD_PROBABILITY
      )

      data.frame(
        fit_id = fit_id_now,
        varying_block = spec_now$varying_block,
        prior_role = spec_now$prior_role,
        parameter = parameter_name,
        n_draws = length(values),
        posterior_mean = mean(values),
        posterior_sd = stats::sd(values),
        HPD80_lower = interval["lower"],
        HPD80_upper = interval["upper"],
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(rbind, rows)
}

parameter_summary <- do.call(
  rbind,
  lapply(sensitivity_results, summarize_parameter_draws)
)
row.names(parameter_summary) <- NULL

write.csv(
  sensitivity_fit_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_fit_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  topk_overlap_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_topk_overlap.csv"
  ),
  row.names = FALSE
)

write.csv(
  parameter_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_parameter_summary.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 13A. EFFECTIVE SAMPLE SIZE DIAGNOSTICS
## ============================================================================

effective_size_matrix <- function(draw_matrix) {
  draw_matrix <- as.matrix(draw_matrix)

  if (
    nrow(draw_matrix) < 2L ||
    ncol(draw_matrix) < 1L ||
    any(!is.finite(draw_matrix))
  ) {
    stop("ESS received an invalid posterior-draw matrix.")
  }

  as.numeric(
    coda::effectiveSize(
      coda::mcmc(draw_matrix)
    )
  )
}

core_parameter_ess_rows <- list()
node_probability_ess_rows <- list()
core_ess_index <- 1L
node_ess_index <- 1L

for (fit_id_now in SENSITIVITY_RUN_ORDER) {
  result_now <- sensitivity_results[[fit_id_now]]
  draws_now <- result_now$fit$draws
  spec_now <- result_now$fit_spec

  beta0_matrix <- matrix(
    as.numeric(draws_now$beta0),
    ncol = 1L,
    dimnames = list(NULL, "beta0")
  )

  eta_matrix <- matrix(
    as.numeric(draws_now$eta),
    ncol = 1L,
    dimnames = list(NULL, "eta")
  )

  core_draw_matrix <- cbind(
    beta0_matrix,
    eta_matrix,
    as.matrix(draws_now$tau),
    as.matrix(draws_now$theta)
  )

  if (
    is.null(colnames(core_draw_matrix)) ||
    anyDuplicated(colnames(core_draw_matrix))
  ) {
    stop(
      "The core posterior-draw columns are not uniquely named for ",
      fit_id_now,
      "."
    )
  }

  core_ess_now <- effective_size_matrix(core_draw_matrix)

  core_parameter_ess_rows[[core_ess_index]] <- data.frame(
    fit_id = fit_id_now,
    varying_block = spec_now$varying_block,
    prior_role = spec_now$prior_role,
    parameter = colnames(core_draw_matrix),
    retained_draws = nrow(core_draw_matrix),
    ESS = core_ess_now,
    ESS_fraction_of_retained =
      core_ess_now / nrow(core_draw_matrix),
    stringsAsFactors = FALSE
  )
  core_ess_index <- core_ess_index + 1L

  latent_draw_matrix <- as.matrix(
    draws_now$latent_probability
  )

  if (!all(dim(latent_draw_matrix) == c(N_RETAINED, N_NODES))) {
    stop(
      "Unexpected latent-probability draw dimensions for ",
      fit_id_now,
      "."
    )
  }

  latent_ess_now <- effective_size_matrix(latent_draw_matrix)

  node_probability_ess_rows[[node_ess_index]] <- data.frame(
    fit_id = fit_id_now,
    varying_block = spec_now$varying_block,
    prior_role = spec_now$prior_role,
    node = seq_len(N_NODES),
    node_group = NODE_GROUP,
    retained_draws = nrow(latent_draw_matrix),
    ESS = latent_ess_now,
    ESS_fraction_of_retained =
      latent_ess_now / nrow(latent_draw_matrix),
    stringsAsFactors = FALSE
  )
  node_ess_index <- node_ess_index + 1L
}

core_parameter_ess <- do.call(
  rbind,
  core_parameter_ess_rows
)
row.names(core_parameter_ess) <- NULL

node_probability_ess <- do.call(
  rbind,
  node_probability_ess_rows
)
row.names(node_probability_ess) <- NULL

node_ess_summary <- do.call(
  rbind,
  lapply(
    split(node_probability_ess, node_probability_ess$fit_id),
    function(df_now) {
      data.frame(
        fit_id = df_now$fit_id[1L],
        varying_block = df_now$varying_block[1L],
        prior_role = df_now$prior_role[1L],
        retained_draws = df_now$retained_draws[1L],
        min_ESS = min(df_now$ESS, na.rm = TRUE),
        q25_ESS = unname(
          stats::quantile(df_now$ESS, 0.25, na.rm = TRUE)
        ),
        median_ESS = stats::median(df_now$ESS, na.rm = TRUE),
        mean_ESS = mean(df_now$ESS, na.rm = TRUE),
        q75_ESS = unname(
          stats::quantile(df_now$ESS, 0.75, na.rm = TRUE)
        ),
        max_ESS = max(df_now$ESS, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  )
)
row.names(node_ess_summary) <- NULL

write.csv(
  core_parameter_ess,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_ESS_core_parameters.csv"
  ),
  row.names = FALSE
)

write.csv(
  node_probability_ess,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_ESS_latent_probabilities_by_node.csv"
  ),
  row.names = FALSE
)

write.csv(
  node_ess_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_ESS_latent_probabilities_summary.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 14. COMMON DEFAULT-BASED ORDERING FOR ALL 15 PANELS
## ============================================================================

baseline_order_reference <- node_summaries_long[
  node_summaries_long$fit_id == "DEFAULT",
  c("node", "node_group", "latent_probability_mean"),
  drop = FALSE
]

if (nrow(baseline_order_reference) != N_NODES) {
  stop("Expected 100 default-fit rows when constructing the plot order.")
}

group_levels_toy <- c(
  "true zero",
  "hidden positive",
  "observed positive"
)

baseline_order_reference$group_order <- match(
  baseline_order_reference$node_group,
  group_levels_toy
)

if (anyNA(baseline_order_reference$group_order)) {
  stop("Unexpected node_group value in the default summaries.")
}

baseline_order_reference <- baseline_order_reference[
  order(
    baseline_order_reference$group_order,
    baseline_order_reference$latent_probability_mean,
    baseline_order_reference$node
  ),
  ,
  drop = FALSE
]

common_plot_order <- baseline_order_reference$node

if (
  length(common_plot_order) != N_NODES ||
  anyDuplicated(common_plot_order) ||
  !setequal(common_plot_order, seq_len(N_NODES))
) {
  stop("The common sensitivity-plot order is invalid.")
}

write.csv(
  data.frame(
    plot_position = seq_len(N_NODES),
    original_node = common_plot_order,
    node_group = baseline_order_reference$node_group,
    default_latent_posterior_mean =
      baseline_order_reference$latent_probability_mean,
    stringsAsFactors = FALSE
  ),
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_common_order.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 15. BUILD THE FINAL 5 X 3 PRIOR-SENSITIVITY FIGURE
## ============================================================================

## This section only rebuilds the figure. It does not rerun the MCMC fits.
if (!exists("OUT_DIR", inherits = TRUE)) OUT_DIR <- file.path(getwd(), "TOY_PU_POE_PRIOR_SENSITIVITY_SEED53")
summary_file <- file.path(OUT_DIR, "TOY_SEED53_PU_POE_PRIOR_SENSITIVITY_COMPLETE_SUMMARIES.RData")
required_objects <- c("node_summaries_long", "sensitivity_prior_bundles", "common_plot_order")
if (any(!vapply(required_objects, exists, logical(1), inherits = TRUE))) {
  if (!file.exists(summary_file)) stop("Completed sensitivity summaries not found: ", summary_file)
  load(summary_file)
}
if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
  stop("Packages 'ggplot2' and 'patchwork' are required.")
}

N_NODES <- length(common_plot_order)
if (N_NODES != 100L || anyDuplicated(common_plot_order) || !setequal(common_plot_order, seq_len(N_NODES))) {
  stop("The common plotting order is invalid.")
}

## --------------------------------------------------------------------------
## Colours
## --------------------------------------------------------------------------

blend_colour <- function(colour, target_colour, amount) {
  amount <- min(max(as.numeric(amount), 0), 1)
  source_rgb <- grDevices::col2rgb(colour)
  target_rgb <- grDevices::col2rgb(target_colour)
  blended_rgb <- (1 - amount) * source_rgb + amount * target_rgb
  grDevices::rgb(blended_rgb[1L, 1L], blended_rgb[2L, 1L], blended_rgb[3L, 1L], maxColorValue = 255)
}

BASE_PUPOE_COLOUR <- "green3"
conservative_base <- blend_colour(BASE_PUPOE_COLOUR, "white", 0.05)
diffused_base <- blend_colour(BASE_PUPOE_COLOUR, "black", 0.05)
column_colours <- c(
  conservative = blend_colour(conservative_base, "yellow", 0.6),
  default = BASE_PUPOE_COLOUR,
  diffuse = blend_colour(diffused_base, "blue", 0.5)
)

## --------------------------------------------------------------------------
## Labels and panel map
## --------------------------------------------------------------------------

baseline_eta_beta <- sensitivity_prior_bundles$DEFAULT$eta_prior$b_eta
baseline_eta_beta_display <- as.numeric(formatC(baseline_eta_beta, format = "f", digits = 2L))
column_order <- c("conservative", "default", "diffuse")
column_headers <- c(
  conservative = "Conservative prior",
  default = "Default prior setting",
  diffuse = "Diffuse prior"
)

panel_specs <- list(
  sigma = list(
    row_label = expression(sigma[j]),
    fit_ids = c("SIGMA_CONSERVATIVE", "DEFAULT", "SIGMA_DIFFUSE"),
    titles = list("Total GP SD = 1", "Total GP SD = 2", "Total GP SD = 4")
  ),
  ell = list(
    row_label = expression(l[j]),
    fit_ids = c("ELL_CONSERVATIVE", "DEFAULT", "ELL_DIFFUSE"),
    titles = list(
      expression("Correlation" %in% "[0.05, 0.25]"),
      expression("Correlation" %in% "[0.10, 0.50]"),
      expression("Correlation" %in% "[0.10, 0.90]")
    )
  ),
  gamma = list(
    row_label = expression(gamma),
    fit_ids = c("GAMMA_CONSERVATIVE", "DEFAULT", "GAMMA_DIFFUSE"),
    titles = list(expression(0.5 * s[gamma]), expression(s[gamma]), expression(2 * s[gamma]))
  ),
  beta0 = list(
    row_label = expression(beta[0]),
    fit_ids = c("BETA0_CONSERVATIVE", "DEFAULT", "BETA0_DIFFUSE"),
    titles = list(
      expression(logit^{-1}*(beta[0]) %in% "[0.005, 0.05]"),
      expression(logit^{-1}*(beta[0]) %in% "[0.01, 0.10]"),
      expression(logit^{-1}*(beta[0]) %in% "[0.02, 0.20]")
    )
  ),
  eta = list(
    row_label = expression(eta),
    fit_ids = c("ETA_CONSERVATIVE", "DEFAULT", "ETA_DIFFUSE"),
    titles = list(
      expression(eta * " ~ " * Beta(1, 5)),
      bquote(eta * " ~ " * Beta(2, .(baseline_eta_beta_display))),
      expression(eta * " ~ " * Beta(1, 1))
    )
  )
)
row_order <- c("sigma", "ell", "gamma", "beta0", "eta")

## Retain a compact map because the original script saves this object later.
plot_grid_map <- do.call(rbind, lapply(row_order, function(row_key_now) {
  row_spec_now <- panel_specs[[row_key_now]]
  data.frame(
    row_key = row_key_now,
    row_label = paste(deparse(row_spec_now$row_label), collapse = ""),
    column_key = column_order,
    fit_id = row_spec_now$fit_ids,
    panel_title = vapply(row_spec_now$titles, function(x) paste(deparse(x), collapse = ""), character(1)),
    stringsAsFactors = FALSE
  )
}))

## --------------------------------------------------------------------------
## Panel function
## --------------------------------------------------------------------------

make_sensitivity_panel <- function(fit_id_now, column_key_now, panel_title_now,
                                   show_x_title = FALSE, show_y_title = FALSE) {
  df <- node_summaries_long[node_summaries_long$fit_id == fit_id_now, , drop = FALSE]
  if (nrow(df) != N_NODES) stop("Expected ", N_NODES, " rows for ", fit_id_now, ".")
  
  df$plot_id <- match(df$node, common_plot_order)
  df <- df[order(df$plot_id), , drop = FALSE]
  df$plot_mean <- df$latent_probability_mean
  df$plot_lower <- df$latent_HPD80_lower
  df$plot_upper <- df$latent_HPD80_upper
  
  if (any(!is.finite(df$plot_mean)) || any(!is.finite(df$plot_lower)) ||
      any(!is.finite(df$plot_upper)) || any(df$plot_lower > df$plot_upper)) {
    stop("Invalid posterior summaries for ", fit_id_now, ".")
  }
  
  colour_now <- unname(column_colours[column_key_now])
  
  ggplot2::ggplot(df, ggplot2::aes(x = plot_id)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = plot_lower, ymax = plot_upper), fill = colour_now, alpha = 0.22) +
    ggplot2::geom_line(ggplot2::aes(y = plot_lower), colour = colour_now, linewidth = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = plot_upper), colour = colour_now, linewidth = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = plot_mean), colour = colour_now, linewidth = 0.65, linetype = 2) +
    ggplot2::geom_vline(xintercept = c(70.5, 80.5), linetype = "dashed", linewidth = 0.33, colour = "grey45") +
    ggplot2::scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100), limits = c(0, 100),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 0.25, 0.50, 0.75, 1.00), limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::labs(
      title = panel_title_now,
      x = if (show_x_title) "Unit, ordered within group" else NULL,
      y = if (show_y_title) "Post. Prob. of T = 1" else NULL
    ) +
    ggplot2::theme_minimal(base_size = 8.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "plain", size = 8.3, hjust = 0.5, margin = ggplot2::margin(b = 2)),
      axis.title = ggplot2::element_text(size = 7.5),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 3)),
      axis.text = ggplot2::element_text(size = 6.8),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey88"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

make_text_strip <- function(label, fontsize) {
  patchwork::wrap_elements(full = grid::textGrob(
    label = label, x = 0.5, y = 0.5,
    gp = grid::gpar(fontface = "bold", fontsize = fontsize)
  ))
}

## --------------------------------------------------------------------------
## Assemble and save
## --------------------------------------------------------------------------

plot_list <- list(
  make_text_strip("", 1),
  make_text_strip(column_headers["conservative"], 12),
  make_text_strip(column_headers["default"], 12),
  make_text_strip(column_headers["diffuse"], 12)
)

for (row_index in seq_along(row_order)) {
  row_spec <- panel_specs[[row_order[row_index]]]
  plot_list[[length(plot_list) + 1L]] <- make_text_strip(row_spec$row_label, 15)
  for (column_index in seq_along(column_order)) {
    plot_list[[length(plot_list) + 1L]] <- make_sensitivity_panel(
      fit_id_now = row_spec$fit_ids[column_index],
      column_key_now = column_order[column_index],
      panel_title_now = row_spec$titles[[column_index]],
      show_x_title = row_index == length(row_order),
      show_y_title = column_index == 1L
    )
  }
}

layout_5x3 <- patchwork::wrap_plots(
  plot_list, ncol = 4, widths = c(0.08, 1, 1, 1),
  heights = c(0.08, 1, 1, 1, 1, 1), byrow = TRUE
) + patchwork::plot_annotation(
  title = "Prior sensitivity analysis for the PU-PoE illustration",
  theme = ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5, margin = ggplot2::margin(b = 4))
  )
)

out_pdf_5x3 <- file.path(OUT_DIR, "PU_PoE_PRIOR_SENSITIVITY_5x3_FINAL_HPD80.pdf")
out_png_5x3 <- file.path(OUT_DIR, "PU_PoE_PRIOR_SENSITIVITY_5x3_FINAL_HPD80.png")
FIGURE_WIDTH_IN <- 8.27
FIGURE_HEIGHT_IN <- 11.69

grDevices::pdf(out_pdf_5x3, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, useDingbats = FALSE)
print(layout_5x3)
grDevices::dev.off()

ggplot2::ggsave(
  filename = out_png_5x3, plot = layout_5x3,
  width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN,
  units = "in", dpi = 300, limitsize = FALSE
)

message("Saved the final prior-sensitivity figure to: ", out_pdf_5x3)

## ============================================================================
## 16. SAVE THE COMPLETE SUMMARY BUNDLE
## ============================================================================

total_runtime_sec <- as.numeric(difftime(
  all_fits_end,
  all_fits_start,
  units = "secs"
))

sensitivity_result_index <- lapply(
  sensitivity_results,
  function(result) {
    list(
      fit_id = result$fit_spec$fit_id,
      fit_spec = result$fit_spec,
      fit_signature = result$fit_signature,
      node_summary = result$node_summary,
      eta_summary = result$eta_summary,
      diagnostics = result$diagnostics,
      prior_calibration = result$prior_calibration,
      completed_at = result$completed_at,
      checkpoint_file = file.path(
        FIT_CHECKPOINT_DIR,
        paste0(result$fit_spec$fit_id, "_complete.rds")
      )
    )
  }
)

run_config <- list(
  script_version = SCRIPT_VERSION,
  created_at = timestamp_now(),
  output_directory = OUT_DIR,
  exact_toy_data_file = TOY_RDS_FILE,
  exact_toy_data_hash = TOY_DATA_HASH,
  toy_seed = TOY_SEED,
  common_mcmc_seed = COMMON_MCMC_SEED,
  n_distinct_fits = length(SENSITIVITY_RUN_ORDER),
  fit_order = SENSITIVITY_RUN_ORDER,
  fit_specs = SENSITIVITY_SPECS,
  MCMC = list(
    blocked_burnin = BLOCKED_BURNIN,
    joint_burnin = JOINT_BURNIN,
    posterior_iterations = POSTERIOR_ITERATIONS,
    thinning = THIN,
    retained_draws = N_RETAINED
  ),
  intervals = list(
    level = HPD_PROBABILITY,
    type = "highest posterior density",
    implementation = "coda::HPDinterval",
    equal_tail_intervals_used = FALSE
  ),
  plotting = list(
    base_colour = BASE_PUPOE_COLOUR,
    conservative_colour = column_colours["conservative"],
    default_colour = column_colours["default"],
    diffuse_colour = column_colours["diffuse"],
    ordering = paste0(
      "Within true-zero, hidden-positive and observed-positive groups, ",
      "order by the default PU-PoE latent posterior mean."
    ),
    figure_pdf = out_pdf_5x3,
    figure_png = out_png_5x3
  ),
  top_k_values = TOP_K_VALUES,
  total_runtime_sec = total_runtime_sec
)

saveRDS(
  run_config,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_run_config.rds"
  )
)

save(
  toy,
  X_cov,
  X_standardization,
  A_network,
  modularity,
  Z_network,
  network_geometry,
  SENSITIVITY_SPECS,
  sensitivity_prior_bundles,
  sensitivity_result_index,
  node_summaries_long,
  eta_summary,
  fit_diagnostics,
  prior_calibration_completed,
  sensitivity_fit_summary,
  topk_overlap_summary,
  parameter_summary,
  core_parameter_ess,
  node_probability_ess,
  node_ess_summary,
  plot_grid_map,
  common_plot_order,
  column_colours,
  run_config,
  file = file.path(
    OUT_DIR,
    "TOY_SEED53_PU_POE_PRIOR_SENSITIVITY_COMPLETE_SUMMARIES.RData"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    OUT_DIR,
    "sessionInfo_toy_seed53_PU_PoE_prior_sensitivity.txt"
  )
)

## ============================================================================
## 17. FINAL CONSOLE OUTPUT
## ============================================================================

cat("\n================ PRIOR SENSITIVITY COMPLETE ================\n")
cat("Toy data seed           :", TOY_SEED, "\n")
cat("Distinct PU-PoE fits    :", length(SENSITIVITY_RUN_ORDER), "\n")
cat("Common MCMC seed        :", COMMON_MCMC_SEED, "\n")
cat("Blocked warm-up/fit     :", BLOCKED_BURNIN, "\n")
cat("Joint warm-up/fit       :", JOINT_BURNIN, "\n")
cat("Posterior iter/fit      :", POSTERIOR_ITERATIONS, "\n")
cat("Thinning                :", THIN, "\n")
cat("Retained draws/fit      :", N_RETAINED, "\n")
cat("Node intervals          : 80% HPD from coda::HPDinterval\n")
cat("Total current runtime   :", format_duration(total_runtime_sec), "\n")

cat("\n================ FIT-LEVEL ROBUSTNESS SUMMARY ================\n")
print(sensitivity_fit_summary, row.names = FALSE)

cat("\n================ ETA POSTERIOR SUMMARY ================\n")
print(eta_summary, row.names = FALSE)

cat("\n================ NODEWISE ESS SUMMARY ================\n")
print(node_ess_summary, row.names = FALSE)

cat("\nFinal 5 x 3 figure:\n")
cat(normalizePath(out_pdf_5x3, winslash = "/", mustWork = FALSE), "\n")
cat("\nAll outputs were written to:\n")
cat(normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE), "\n")
cat("\n================ DONE ================\n")





###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################









## ============================================================================
## 18. JOINT PRIOR STRESS TEST: TWO EXTRA FITS + COMBINED 6 X 3 FIGURE
## ----------------------------------------------------------------------------
## Append this patch after the current Section 17.
##
## When the full script is sourced again, the original 11 fits are loaded from
## their checkpoints and only the following two new fits are computed:
##
##   JOINT_CONSERVATIVE : all five conservative prior specifications together
##   JOINT_DIFFUSE      : all five diffuse prior specifications together
##
## The patch then rebuilds every numerical summary for all 13 fits and creates
## a 6 x 3 figure. The sixth row is labelled "Joint priors" vertically.
## ============================================================================

required_joint_objects <- c(
  "SENSITIVITY_SPECS",
  "SENSITIVITY_RUN_ORDER",
  "sensitivity_prior_bundles",
  "sensitivity_results",
  "sensitivity_fit_summary",
  "topk_overlap_summary",
  "parameter_summary",
  "core_parameter_ess",
  "node_probability_ess",
  "node_ess_summary",
  "node_summaries_long",
  "eta_summary",
  "fit_diagnostics",
  "prior_calibration_completed",
  "common_plot_order",
  "panel_specs",
  "row_order",
  "column_order",
  "column_headers",
  "column_colours",
  "make_sensitivity_spec",
  "calibrate_sensitivity_fit_priors",
  "run_or_load_sensitivity_fit",
  "summarize_parameter_draws",
  "effective_size_matrix",
  "ordered_top_nodes"
)

missing_joint_objects <- required_joint_objects[
  !vapply(required_joint_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_joint_objects) > 0L) {
  stop(
    "Run the original sensitivity script before this patch. Missing object(s): ",
    paste(missing_joint_objects, collapse = ", ")
  )
}

## ============================================================================
## 18A. DEFINE AND CALIBRATE THE TWO JOINT PRIOR SPECIFICATIONS
## ============================================================================

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

for (fit_id_now in names(JOINT_STRESS_SPECS)) {
  SENSITIVITY_SPECS[[fit_id_now]] <- JOINT_STRESS_SPECS[[fit_id_now]]
  sensitivity_prior_bundles[[fit_id_now]] <-
    calibrate_sensitivity_fit_priors(
      fit_spec = JOINT_STRESS_SPECS[[fit_id_now]],
      X_cov = X_cov,
      Z_network = Z_network,
      network_geometry = network_geometry
    )
}

original_run_order <- setdiff(
  SENSITIVITY_RUN_ORDER,
  names(JOINT_STRESS_SPECS)
)

SENSITIVITY_RUN_ORDER_ALL <- c(
  original_run_order,
  names(JOINT_STRESS_SPECS)
)

## ============================================================================
## 18B. RUN OR LOAD ONLY THE TWO EXTRA FITS
## ============================================================================

joint_fit_start <- Sys.time()

joint_sensitivity_results <- vector(
  "list",
  length(JOINT_STRESS_SPECS)
)
names(joint_sensitivity_results) <- names(JOINT_STRESS_SPECS)

for (fit_id_now in names(JOINT_STRESS_SPECS)) {
  joint_sensitivity_results[[fit_id_now]] <-
    run_or_load_sensitivity_fit(
      JOINT_STRESS_SPECS[[fit_id_now]]
    )
  invisible(gc())
}

joint_fit_end <- Sys.time()

if (any(!vapply(
  joint_sensitivity_results,
  function(x) is.list(x) && isTRUE(x$complete),
  logical(1)
))) {
  stop("At least one joint prior stress-test fit did not complete.")
}

sensitivity_results_all <- sensitivity_results

for (fit_id_now in names(joint_sensitivity_results)) {
  sensitivity_results_all[[fit_id_now]] <-
    joint_sensitivity_results[[fit_id_now]]
}

sensitivity_results_all <- sensitivity_results_all[
  SENSITIVITY_RUN_ORDER_ALL
]

## ============================================================================
## 18C. COMBINE THE EXISTING 11 FITS WITH THE TWO JOINT FITS
## ============================================================================

joint_node_summaries <- do.call(
  rbind,
  lapply(joint_sensitivity_results, `[[`, "node_summary")
)
row.names(joint_node_summaries) <- NULL

joint_eta_summary_rows <- do.call(
  rbind,
  lapply(joint_sensitivity_results, `[[`, "eta_summary")
)
row.names(joint_eta_summary_rows) <- NULL

joint_fit_diagnostics <- do.call(
  rbind,
  lapply(joint_sensitivity_results, `[[`, "diagnostics")
)
row.names(joint_fit_diagnostics) <- NULL

joint_prior_calibration <- do.call(
  rbind,
  lapply(joint_sensitivity_results, `[[`, "prior_calibration")
)
row.names(joint_prior_calibration) <- NULL

node_summaries_long_all <- rbind(
  node_summaries_long,
  joint_node_summaries
)

eta_summary_all <- rbind(
  eta_summary,
  joint_eta_summary_rows
)

fit_diagnostics_all <- rbind(
  fit_diagnostics,
  joint_fit_diagnostics
)

prior_calibration_completed_all <- rbind(
  prior_calibration_completed,
  joint_prior_calibration
)

fit_order_index <- match(
  node_summaries_long_all$fit_id,
  SENSITIVITY_RUN_ORDER_ALL
)

node_summaries_long_all <- node_summaries_long_all[
  order(fit_order_index, node_summaries_long_all$node),
  ,
  drop = FALSE
]
row.names(node_summaries_long_all) <- NULL

eta_summary_all <- eta_summary_all[
  order(match(eta_summary_all$fit_id, SENSITIVITY_RUN_ORDER_ALL)),
  ,
  drop = FALSE
]
row.names(eta_summary_all) <- NULL

fit_diagnostics_all <- fit_diagnostics_all[
  order(match(fit_diagnostics_all$fit_id, SENSITIVITY_RUN_ORDER_ALL)),
  ,
  drop = FALSE
]
row.names(fit_diagnostics_all) <- NULL

prior_calibration_completed_all <- prior_calibration_completed_all[
  order(
    match(
      prior_calibration_completed_all$fit_id,
      SENSITIVITY_RUN_ORDER_ALL
    ),
    prior_calibration_completed_all$expert
  ),
  ,
  drop = FALSE
]
row.names(prior_calibration_completed_all) <- NULL

## ============================================================================
## 18D. FIT-LEVEL ROBUSTNESS AND TOP-K SUMMARIES FOR THE TWO EXTRA FITS
## ============================================================================

baseline_nodes_all <- node_summaries_long_all[
  node_summaries_long_all$fit_id == "DEFAULT",
  c("node", "latent_probability_mean"),
  drop = FALSE
]

baseline_nodes_all <- baseline_nodes_all[
  order(baseline_nodes_all$node),
  ,
  drop = FALSE
]

if (nrow(baseline_nodes_all) != N_NODES) {
  stop("The default fit did not return exactly 100 node summaries.")
}

baseline_probability_all <- baseline_nodes_all$latent_probability_mean
names(baseline_probability_all) <- baseline_nodes_all$node

joint_fit_summary_rows <- list()
joint_topk_rows <- list()
joint_fit_summary_index <- 1L
joint_topk_index <- 1L

for (fit_id_now in names(JOINT_STRESS_SPECS)) {
  result_now <- sensitivity_results_all[[fit_id_now]]
  spec_now <- result_now$fit_spec
  prior_now <- result_now$prior_bundle
  
  node_now <- result_now$node_summary[
    order(result_now$node_summary$node),
    ,
    drop = FALSE
  ]
  
  probability_now <- node_now$latent_probability_mean
  
  spearman_now <- suppressWarnings(
    stats::cor(
      probability_now,
      baseline_probability_all,
      method = "spearman"
    )
  )
  
  eta_now <- result_now$eta_summary[1L, , drop = FALSE]
  diagnostics_now <- result_now$diagnostics[1L, , drop = FALSE]
  
  joint_fit_summary_rows[[joint_fit_summary_index]] <- data.frame(
    fit_id = fit_id_now,
    varying_block = spec_now$varying_block,
    prior_role = spec_now$prior_role,
    total_sd_center = spec_now$total_sd_center,
    corr_pct_low = spec_now$corr_pct_range[1L],
    corr_pct_high = spec_now$corr_pct_range[2L],
    mean_rms_95_mult = spec_now$mean_rms_95_mult,
    beta0_prev_low = spec_now$beta0_prev_range[1L],
    beta0_prev_high = spec_now$beta0_prev_range[2L],
    beta0_prior_mean = prior_now$beta0_prior$mean,
    beta0_prior_sd = prior_now$beta0_prior$sd,
    eta_prior_alpha = prior_now$eta_prior$a_eta,
    eta_prior_beta = prior_now$eta_prior$b_eta,
    eta_posterior_mean = eta_now$eta_posterior_mean,
    eta_posterior_sd = eta_now$eta_posterior_sd,
    eta_HPD80_lower = eta_now$eta_HPD80_lower,
    eta_HPD80_upper = eta_now$eta_HPD80_upper,
    spearman_with_default = spearman_now,
    mean_absolute_probability_difference = mean(
      abs(probability_now - baseline_probability_all)
    ),
    max_absolute_probability_difference = max(
      abs(probability_now - baseline_probability_all)
    ),
    mean_HPD80_width = mean(node_now$latent_HPD80_width),
    overall_AUC = diagnostics_now$overall_AUC_latent_mean,
    hidden_AUC = diagnostics_now$hidden_AUC_latent_mean,
    runtime_sec = result_now$fit$runtime_sec,
    stringsAsFactors = FALSE
  )
  
  joint_fit_summary_index <- joint_fit_summary_index + 1L
  
  for (k_now in TOP_K_VALUES) {
    current_top <- ordered_top_nodes(
      probability_now,
      k_now
    )
    
    baseline_top <- ordered_top_nodes(
      baseline_probability_all,
      k_now
    )
    
    overlap_count <- length(
      intersect(current_top, baseline_top)
    )
    
    joint_topk_rows[[joint_topk_index]] <- data.frame(
      fit_id = fit_id_now,
      varying_block = spec_now$varying_block,
      prior_role = spec_now$prior_role,
      k = k_now,
      overlap_count = overlap_count,
      overlap_fraction = overlap_count / k_now,
      stringsAsFactors = FALSE
    )
    
    joint_topk_index <- joint_topk_index + 1L
  }
}

joint_fit_summary <- do.call(
  rbind,
  joint_fit_summary_rows
)
row.names(joint_fit_summary) <- NULL

joint_topk_summary <- do.call(
  rbind,
  joint_topk_rows
)
row.names(joint_topk_summary) <- NULL

sensitivity_fit_summary_all <- rbind(
  sensitivity_fit_summary,
  joint_fit_summary
)

sensitivity_fit_summary_all <- sensitivity_fit_summary_all[
  order(
    match(
      sensitivity_fit_summary_all$fit_id,
      SENSITIVITY_RUN_ORDER_ALL
    )
  ),
  ,
  drop = FALSE
]
row.names(sensitivity_fit_summary_all) <- NULL

topk_overlap_summary_all <- rbind(
  topk_overlap_summary,
  joint_topk_summary
)

topk_overlap_summary_all <- topk_overlap_summary_all[
  order(
    match(
      topk_overlap_summary_all$fit_id,
      SENSITIVITY_RUN_ORDER_ALL
    ),
    topk_overlap_summary_all$k
  ),
  ,
  drop = FALSE
]
row.names(topk_overlap_summary_all) <- NULL

## ============================================================================
## 18E. PARAMETER POSTERIORS AND ESS FOR THE TWO EXTRA FITS
## ============================================================================

joint_parameter_summary <- do.call(
  rbind,
  lapply(
    joint_sensitivity_results,
    summarize_parameter_draws
  )
)
row.names(joint_parameter_summary) <- NULL

parameter_summary_all <- rbind(
  parameter_summary,
  joint_parameter_summary
)

parameter_summary_all <- parameter_summary_all[
  order(
    match(
      parameter_summary_all$fit_id,
      SENSITIVITY_RUN_ORDER_ALL
    ),
    parameter_summary_all$parameter
  ),
  ,
  drop = FALSE
]
row.names(parameter_summary_all) <- NULL

joint_core_parameter_ess_rows <- list()
joint_node_probability_ess_rows <- list()
joint_core_index <- 1L
joint_node_index <- 1L

for (fit_id_now in names(JOINT_STRESS_SPECS)) {
  result_now <- sensitivity_results_all[[fit_id_now]]
  draws_now <- result_now$fit$draws
  spec_now <- result_now$fit_spec
  
  beta0_matrix <- matrix(
    as.numeric(draws_now$beta0),
    ncol = 1L,
    dimnames = list(NULL, "beta0")
  )
  
  eta_matrix <- matrix(
    as.numeric(draws_now$eta),
    ncol = 1L,
    dimnames = list(NULL, "eta")
  )
  
  core_draw_matrix <- cbind(
    beta0_matrix,
    eta_matrix,
    as.matrix(draws_now$tau),
    as.matrix(draws_now$theta)
  )
  
  if (
    is.null(colnames(core_draw_matrix)) ||
    anyDuplicated(colnames(core_draw_matrix))
  ) {
    stop(
      "The core posterior-draw columns are not uniquely named for ",
      fit_id_now,
      "."
    )
  }
  
  core_ess_now <- effective_size_matrix(
    core_draw_matrix
  )
  
  joint_core_parameter_ess_rows[[joint_core_index]] <- data.frame(
    fit_id = fit_id_now,
    varying_block = spec_now$varying_block,
    prior_role = spec_now$prior_role,
    parameter = colnames(core_draw_matrix),
    retained_draws = nrow(core_draw_matrix),
    ESS = core_ess_now,
    ESS_fraction_of_retained =
      core_ess_now / nrow(core_draw_matrix),
    stringsAsFactors = FALSE
  )
  
  joint_core_index <- joint_core_index + 1L
  
  latent_draw_matrix <- as.matrix(
    draws_now$latent_probability
  )
  
  if (!all(dim(latent_draw_matrix) == c(N_RETAINED, N_NODES))) {
    stop(
      "Unexpected latent-probability draw dimensions for ",
      fit_id_now,
      "."
    )
  }
  
  latent_ess_now <- effective_size_matrix(
    latent_draw_matrix
  )
  
  joint_node_probability_ess_rows[[joint_node_index]] <- data.frame(
    fit_id = fit_id_now,
    varying_block = spec_now$varying_block,
    prior_role = spec_now$prior_role,
    node = seq_len(N_NODES),
    node_group = NODE_GROUP,
    retained_draws = nrow(latent_draw_matrix),
    ESS = latent_ess_now,
    ESS_fraction_of_retained =
      latent_ess_now / nrow(latent_draw_matrix),
    stringsAsFactors = FALSE
  )
  
  joint_node_index <- joint_node_index + 1L
}

joint_core_parameter_ess <- do.call(
  rbind,
  joint_core_parameter_ess_rows
)
row.names(joint_core_parameter_ess) <- NULL

joint_node_probability_ess <- do.call(
  rbind,
  joint_node_probability_ess_rows
)
row.names(joint_node_probability_ess) <- NULL

joint_node_ess_summary <- do.call(
  rbind,
  lapply(
    split(
      joint_node_probability_ess,
      joint_node_probability_ess$fit_id
    ),
    function(df_now) {
      data.frame(
        fit_id = df_now$fit_id[1L],
        varying_block = df_now$varying_block[1L],
        prior_role = df_now$prior_role[1L],
        retained_draws = df_now$retained_draws[1L],
        min_ESS = min(df_now$ESS, na.rm = TRUE),
        q25_ESS = unname(
          stats::quantile(
            df_now$ESS,
            0.25,
            na.rm = TRUE
          )
        ),
        median_ESS = stats::median(
          df_now$ESS,
          na.rm = TRUE
        ),
        mean_ESS = mean(
          df_now$ESS,
          na.rm = TRUE
        ),
        q75_ESS = unname(
          stats::quantile(
            df_now$ESS,
            0.75,
            na.rm = TRUE
          )
        ),
        max_ESS = max(
          df_now$ESS,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)
row.names(joint_node_ess_summary) <- NULL

core_parameter_ess_all <- rbind(
  core_parameter_ess,
  joint_core_parameter_ess
)

core_parameter_ess_all <- core_parameter_ess_all[
  order(
    match(
      core_parameter_ess_all$fit_id,
      SENSITIVITY_RUN_ORDER_ALL
    ),
    core_parameter_ess_all$parameter
  ),
  ,
  drop = FALSE
]
row.names(core_parameter_ess_all) <- NULL

node_probability_ess_all <- rbind(
  node_probability_ess,
  joint_node_probability_ess
)

node_probability_ess_all <- node_probability_ess_all[
  order(
    match(
      node_probability_ess_all$fit_id,
      SENSITIVITY_RUN_ORDER_ALL
    ),
    node_probability_ess_all$node
  ),
  ,
  drop = FALSE
]
row.names(node_probability_ess_all) <- NULL

node_ess_summary_all <- rbind(
  node_ess_summary,
  joint_node_ess_summary
)

node_ess_summary_all <- node_ess_summary_all[
  order(
    match(
      node_ess_summary_all$fit_id,
      SENSITIVITY_RUN_ORDER_ALL
    )
  ),
  ,
  drop = FALSE
]
row.names(node_ess_summary_all) <- NULL

## ============================================================================
## 18F. COMPACT THREE-WAY SUMMARY FOR THE JOINT STRESS TEST
## ============================================================================

JOINT_COMPARISON_IDS <- c(
  "JOINT_CONSERVATIVE",
  "DEFAULT",
  "JOINT_DIFFUSE"
)

joint_eta_comparison <- do.call(
  rbind,
  lapply(
    JOINT_COMPARISON_IDS,
    function(fit_id_now) {
      eta_row <- eta_summary_all[
        eta_summary_all$fit_id == fit_id_now,
        ,
        drop = FALSE
      ]
      
      if (nrow(eta_row) != 1L) {
        stop(
          "Expected exactly one eta-summary row for ",
          fit_id_now,
          "."
        )
      }
      
      result_now <- sensitivity_results_all[[fit_id_now]]
      eta_prior_now <- result_now$prior_bundle$eta_prior
      
      display_role <- switch(
        fit_id_now,
        JOINT_CONSERVATIVE = "All conservative",
        DEFAULT = "Default",
        JOINT_DIFFUSE = "All diffuse"
      )
      
      data.frame(
        fit_id = fit_id_now,
        specification = display_role,
        eta_prior_alpha = eta_prior_now$a_eta,
        eta_prior_beta = eta_prior_now$b_eta,
        eta_prior_mean =
          eta_prior_now$a_eta /
          (eta_prior_now$a_eta + eta_prior_now$b_eta),
        eta_true = eta_row$eta_true,
        eta_posterior_mean = eta_row$eta_posterior_mean,
        eta_posterior_sd = eta_row$eta_posterior_sd,
        eta_HPD80_lower = eta_row$eta_HPD80_lower,
        eta_HPD80_upper = eta_row$eta_HPD80_upper,
        eta_HPD80_width = eta_row$eta_HPD80_width,
        stringsAsFactors = FALSE
      )
    }
  )
)
row.names(joint_eta_comparison) <- NULL

joint_stress_fit_comparison <- sensitivity_fit_summary_all[
  match(
    JOINT_COMPARISON_IDS,
    sensitivity_fit_summary_all$fit_id
  ),
  ,
  drop = FALSE
]
row.names(joint_stress_fit_comparison) <- NULL

## ============================================================================
## 18G. WRITE ALL 13-FIT COMPARISON TABLES
## ============================================================================

write.csv(
  node_summaries_long_all,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_node_summaries_HPD80.csv"
  ),
  row.names = FALSE
)

write.csv(
  eta_summary_all,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_eta_HPD80.csv"
  ),
  row.names = FALSE
)

write.csv(
  fit_diagnostics_all,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_MCMC_diagnostics.csv"
  ),
  row.names = FALSE
)

write.csv(
  prior_calibration_completed_all,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_active_priors.csv"
  ),
  row.names = FALSE
)

write.csv(
  sensitivity_fit_summary_all,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_fit_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  topk_overlap_summary_all,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_topk_overlap.csv"
  ),
  row.names = FALSE
)

write.csv(
  parameter_summary_all,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_parameter_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  core_parameter_ess_all,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_ESS_core_parameters.csv"
  ),
  row.names = FALSE
)

write.csv(
  node_probability_ess_all,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_ESS_latent_probabilities_by_node.csv"
  ),
  row.names = FALSE
)

write.csv(
  node_ess_summary_all,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_ESS_latent_probabilities_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  joint_eta_comparison,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_JOINT_PRIOR_STRESS_TEST_eta_HPD80.csv"
  ),
  row.names = FALSE
)

write.csv(
  joint_stress_fit_comparison,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_JOINT_PRIOR_STRESS_TEST_fit_comparison.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 18H. BUILD THE FINAL 6 X 3 FIGURE WITH THE SAME CURRENT PARAMETERS
## ============================================================================

panel_specs_6x3 <- panel_specs

panel_specs_6x3$joint <- list(
  row_label = "Joint priors",
  fit_ids = c(
    "JOINT_CONSERVATIVE",
    "DEFAULT",
    "JOINT_DIFFUSE"
  ),
  titles = list(
    "All conservative",
    "Default",
    "All diffuse"
  )
)

row_order_6x3 <- c(
  row_order,
  "joint"
)

plot_grid_map_6x3 <- do.call(
  rbind,
  lapply(
    row_order_6x3,
    function(row_key_now) {
      row_spec_now <- panel_specs_6x3[[row_key_now]]
      
      data.frame(
        row_key = row_key_now,
        row_label = paste(
          deparse(row_spec_now$row_label),
          collapse = ""
        ),
        column_key = column_order,
        fit_id = row_spec_now$fit_ids,
        panel_title = vapply(
          row_spec_now$titles,
          function(x) paste(deparse(x), collapse = ""),
          character(1)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  plot_grid_map_6x3,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_6x3_plot_grid_map.csv"
  ),
  row.names = FALSE
)

make_sensitivity_panel_6x3 <- function(
    fit_id_now,
    column_key_now,
    panel_title_now,
    show_x_title = FALSE,
    show_y_title = FALSE
) {
  df <- node_summaries_long_all[
    node_summaries_long_all$fit_id == fit_id_now,
    ,
    drop = FALSE
  ]
  
  if (nrow(df) != N_NODES) {
    stop(
      "Expected ",
      N_NODES,
      " rows for ",
      fit_id_now,
      "."
    )
  }
  
  df$plot_id <- match(
    df$node,
    common_plot_order
  )
  
  df <- df[
    order(df$plot_id),
    ,
    drop = FALSE
  ]
  
  df$plot_mean <- df$latent_probability_mean
  df$plot_lower <- df$latent_HPD80_lower
  df$plot_upper <- df$latent_HPD80_upper
  
  if (
    any(!is.finite(df$plot_mean)) ||
    any(!is.finite(df$plot_lower)) ||
    any(!is.finite(df$plot_upper)) ||
    any(df$plot_lower > df$plot_upper)
  ) {
    stop(
      "Invalid posterior summaries for ",
      fit_id_now,
      "."
    )
  }
  
  colour_now <- unname(
    column_colours[column_key_now]
  )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = plot_id)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = plot_lower,
        ymax = plot_upper
      ),
      fill = colour_now,
      alpha = 0.22
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = plot_lower),
      colour = colour_now,
      linewidth = 0.35
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = plot_upper),
      colour = colour_now,
      linewidth = 0.35
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = plot_mean),
      colour = colour_now,
      linewidth = 0.65,
      linetype = 2
    ) +
    ggplot2::geom_vline(
      xintercept = c(70.5, 80.5),
      linetype = "dashed",
      linewidth = 0.33,
      colour = "grey45"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      limits = c(0, 100),
      expand = ggplot2::expansion(
        mult = c(0, 0)
      )
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      limits = c(0, 1),
      expand = ggplot2::expansion(
        mult = c(0.01, 0.02)
      )
    ) +
    ggplot2::labs(
      title = panel_title_now,
      x = if (show_x_title) {
        "Unit, ordered within group"
      } else {
        NULL
      },
      y = if (show_y_title) {
        "Post. Prob. of T = 1"
      } else {
        NULL
      }
    ) +
    ggplot2::theme_minimal(
      base_size = 8.5
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "plain",
        size = 8.3,
        hjust = 0.5,
        margin = ggplot2::margin(b = 2)
      ),
      axis.title = ggplot2::element_text(
        size = 7.5
      ),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 3)
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 3)
      ),
      axis.text = ggplot2::element_text(
        size = 6.8
      ),
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.25,
        colour = "grey88"
      ),
      plot.margin = ggplot2::margin(
        4,
        4,
        4,
        4
      )
    )
}

make_vertical_text_strip <- function(
    label,
    fontsize = 11
) {
  patchwork::wrap_elements(
    full = grid::textGrob(
      label = label,
      x = 0.5,
      y = 0.5,
      rot = 90,
      gp = grid::gpar(
        fontface = "bold",
        fontsize = fontsize
      )
    )
  )
}

plot_list_6x3 <- list(
  make_text_strip("", 1),
  make_text_strip(column_headers["conservative"], 12),
  make_text_strip(column_headers["default"], 12),
  make_text_strip(column_headers["diffuse"], 12)
)

for (row_index in seq_along(row_order_6x3)) {
  row_key_now <- row_order_6x3[row_index]
  row_spec_now <- panel_specs_6x3[[row_key_now]]
  
  if (identical(row_key_now, "joint")) {
    plot_list_6x3[[length(plot_list_6x3) + 1L]] <-
      make_vertical_text_strip(
        "Joint priors",
        fontsize = 11
      )
  } else {
    plot_list_6x3[[length(plot_list_6x3) + 1L]] <-
      make_text_strip(
        row_spec_now$row_label,
        15
      )
  }
  
  for (column_index in seq_along(column_order)) {
    plot_list_6x3[[length(plot_list_6x3) + 1L]] <-
      make_sensitivity_panel_6x3(
        fit_id_now =
          row_spec_now$fit_ids[column_index],
        column_key_now =
          column_order[column_index],
        panel_title_now =
          row_spec_now$titles[[column_index]],
        show_x_title =
          row_index == length(row_order_6x3),
        show_y_title =
          column_index == 1L
      )
  }
}

layout_6x3 <- patchwork::wrap_plots(
  plot_list_6x3,
  ncol = 4,
  widths = c(0.08, 1, 1, 1),
  heights = c(0.08, 1, 1, 1, 1, 1, 1),
  byrow = TRUE
) +
  patchwork::plot_annotation(
    title = "Prior sensitivity analysis for the PU-PoE illustration",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14,
        hjust = 0.5,
        margin = ggplot2::margin(b = 4)
      )
    )
  )

out_pdf_6x3 <- file.path(
  OUT_DIR,
  "PU_PoE_PRIOR_SENSITIVITY_6x3_JOINT_STRESS_HPD80.pdf"
)

out_png_6x3 <- file.path(
  OUT_DIR,
  "PU_PoE_PRIOR_SENSITIVITY_6x3_JOINT_STRESS_HPD80.png"
)

FIGURE_WIDTH_IN_6x3 <- 8.27
FIGURE_HEIGHT_IN_6x3 <- 11.69

grDevices::pdf(
  out_pdf_6x3,
  width = FIGURE_WIDTH_IN_6x3,
  height = FIGURE_HEIGHT_IN_6x3,
  useDingbats = FALSE
)

print(layout_6x3)

grDevices::dev.off()

ggplot2::ggsave(
  filename = out_png_6x3,
  plot = layout_6x3,
  width = FIGURE_WIDTH_IN_6x3,
  height = FIGURE_HEIGHT_IN_6x3,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

## ============================================================================
## 18I. SAVE THE COMBINED 13-FIT SUMMARY BUNDLE
## ============================================================================

sensitivity_result_index_all <- lapply(
  sensitivity_results_all,
  function(result) {
    list(
      fit_id = result$fit_spec$fit_id,
      fit_spec = result$fit_spec,
      fit_signature = result$fit_signature,
      node_summary = result$node_summary,
      eta_summary = result$eta_summary,
      diagnostics = result$diagnostics,
      prior_calibration = result$prior_calibration,
      completed_at = result$completed_at,
      checkpoint_file = file.path(
        FIT_CHECKPOINT_DIR,
        paste0(
          result$fit_spec$fit_id,
          "_complete.rds"
        )
      )
    )
  }
)

joint_patch_runtime_sec <- as.numeric(
  difftime(
    joint_fit_end,
    joint_fit_start,
    units = "secs"
  )
)

joint_stress_run_config <- list(
  script_version = SCRIPT_VERSION,
  created_at = timestamp_now(),
  output_directory = OUT_DIR,
  common_mcmc_seed = COMMON_MCMC_SEED,
  n_distinct_fits = length(SENSITIVITY_RUN_ORDER_ALL),
  fit_order = SENSITIVITY_RUN_ORDER_ALL,
  joint_specs = JOINT_STRESS_SPECS,
  joint_patch_runtime_sec = joint_patch_runtime_sec,
  MCMC = list(
    blocked_burnin = BLOCKED_BURNIN,
    joint_burnin = JOINT_BURNIN,
    posterior_iterations = POSTERIOR_ITERATIONS,
    thinning = THIN,
    retained_draws = N_RETAINED
  ),
  intervals = list(
    level = HPD_PROBABILITY,
    type = "highest posterior density",
    implementation = "coda::HPDinterval"
  ),
  plotting = list(
    row_label = "Joint priors",
    panel_titles = c(
      "All conservative",
      "Default",
      "All diffuse"
    ),
    figure_pdf = out_pdf_6x3,
    figure_png = out_png_6x3
  )
)

saveRDS(
  joint_stress_run_config,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_13fits_joint_stress_run_config.rds"
  )
)

save(
  JOINT_STRESS_SPECS,
  SENSITIVITY_RUN_ORDER_ALL,
  sensitivity_prior_bundles,
  sensitivity_result_index_all,
  node_summaries_long_all,
  eta_summary_all,
  fit_diagnostics_all,
  prior_calibration_completed_all,
  sensitivity_fit_summary_all,
  topk_overlap_summary_all,
  parameter_summary_all,
  core_parameter_ess_all,
  node_probability_ess_all,
  node_ess_summary_all,
  joint_eta_comparison,
  joint_stress_fit_comparison,
  plot_grid_map_6x3,
  common_plot_order,
  column_colours,
  joint_stress_run_config,
  file = file.path(
    OUT_DIR,
    "TOY_SEED53_PU_POE_PRIOR_SENSITIVITY_13FITS_JOINT_STRESS_SUMMARIES.RData"
  )
)

## ============================================================================
## 18J. FINAL COMBINED CONSOLE OUTPUT
## ============================================================================

cat(
  "\n================ 13-FIT PRIOR SENSITIVITY COMPLETE ================\n"
)
cat(
  "Original one-at-a-time fits :",
  length(original_run_order),
  "\n"
)
cat(
  "Extra joint fits            :",
  length(JOINT_STRESS_SPECS),
  "\n"
)
cat(
  "Total distinct fits         :",
  length(SENSITIVITY_RUN_ORDER_ALL),
  "\n"
)
cat(
  "Joint-fit current runtime   :",
  format_duration(joint_patch_runtime_sec),
  "\n"
)

cat(
  "\n================ ALL 13 FIT-LEVEL ROBUSTNESS RESULTS ================\n"
)
print(
  sensitivity_fit_summary_all,
  row.names = FALSE
)

cat(
  "\n================ JOINT STRESS-TEST COMPARISON ================\n"
)
print(
  joint_stress_fit_comparison,
  row.names = FALSE
)

cat(
  "\n================ ALL 13 ETA POSTERIOR RESULTS ================\n"
)
print(
  eta_summary_all,
  row.names = FALSE
)

cat(
  "\n================ JOINT ETA COMPARISON ================\n"
)
print(
  joint_eta_comparison,
  row.names = FALSE
)

cat(
  "\n================ ALL 13 MCMC DIAGNOSTICS ================\n"
)
print(
  fit_diagnostics_all,
  row.names = FALSE
)

cat(
  "\n================ ALL 13 PARAMETER POSTERIOR SUMMARIES ================\n"
)
print(
  parameter_summary_all,
  row.names = FALSE
)

cat(
  "\n================ ALL 13 CORE-PARAMETER ESS RESULTS ================\n"
)
print(
  core_parameter_ess_all,
  row.names = FALSE
)

cat(
  "\n================ ALL 13 NODEWISE ESS SUMMARIES ================\n"
)
print(
  node_ess_summary_all,
  row.names = FALSE
)

cat(
  "\n================ ALL 13 TOP-K OVERLAPS ================\n"
)
print(
  topk_overlap_summary_all,
  row.names = FALSE
)

cat("\nFinal 6 x 3 figure:\n")
cat(
  normalizePath(
    out_pdf_6x3,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)

cat("\nCombined 13-fit summaries:\n")
cat(
  normalizePath(
    file.path(
      OUT_DIR,
      "TOY_SEED53_PU_POE_PRIOR_SENSITIVITY_13FITS_JOINT_STRESS_SUMMARIES.RData"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)

cat(
  "\n================ JOINT PRIOR STRESS TEST DONE ================\n"
)








