## Extracted from GP_four_model_illustrtation.R. See docs/SOURCE_MAP.csv.
## Original mathematical operations and their order are retained.

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

