## Run from project root: Rscript --vanilla tests/test_inputs.R
## Checks the synthetic data, the prior calibration and the likelihood.
stopifnot(file.exists(".here"))
wf <- new.env(parent = as.environment("package:stats"))
sys.source("R/illustration/workflow.R", envir = wf)
wf$check_illustration_dependencies()

e <- wf$new_illustration_environment(
  getwd(),
  "quick",
  tempfile("illustration_test_")
)

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
a <- e$generate_exact_toy_seed53()
b <- e$generate_exact_toy_seed53()

stopifnot(
  identical(a, b),
  sum(a$T) == 30L,
  sum(a$Y) == 20L,
  sum(a$T == 1L & a$Y == 0L) == 10L,
  identical(a$A, t(a$A)),
  all(diag(a$A) == 0L)
)

X <- e$standardize_columns(a$X)$scaled
m1 <- e$modularity_features_single_network(a$A, e$MODULARITY_SEED)
m2 <- e$modularity_features_single_network(a$A, e$MODULARITY_SEED)
stopifnot(identical(m1, m2), max(abs(colMeans(X))) < 1e-12)
g <- e$prepare_network_components_for_blocks(a$A)

p1 <- lapply(
  e$MODEL_SPECS,
  e$calibrate_model_priors,
  X_cov = X,
  Z_network = m1$Z,
  network_geometry = g
)

p2 <- lapply(
  e$MODEL_SPECS,
  e$calibrate_model_priors,
  X_cov = X,
  Z_network = m1$Z,
  network_geometry = g
)

stopifnot(identical(p1, p2))

for (p in p1) {
  stopifnot(
    all(is.finite(p$tau_scale)),
    all(p$tau_scale > 0),
    abs(sum(p$effective_sd_center^2) - e$TOTAL_SD_CENTER^2) < 1e-12
  )
}

f <- c(-4, -1, 0, 1, 4)
y <- c(0, 0, 1, 1, 0)

stopifnot(
  abs(e$pu_loglik_vec(f, y, 0) - e$bernoulli_loglik_vec(f, y)) < 1e-10,
  max(abs(e$pu_gradloglik_vec(f, y, 0) - e$bernoulli_gradloglik_vec(f, y))) <
    1e-10
)

eta <- 0.2
h <- 1e-5

fd <- vapply(
  seq_along(f),
  function(i) {
    plus <- minus <- f
    plus[i] <- plus[i] + h
    minus[i] <- minus[i] - h
    (e$pu_loglik_vec(plus, y, eta) - e$pu_loglik_vec(minus, y, eta)) / (2 * h)
  },
  numeric(1)
)

stopifnot(max(abs(fd - e$pu_gradloglik_vec(f, y, eta))) < 1e-7)

cat(
  "PASS: data, modularity features, prior calibration, variance budget and likelihood.\n"
)
