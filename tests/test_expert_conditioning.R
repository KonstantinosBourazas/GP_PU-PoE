## Run from the open project: source("tests/test_expert_conditioning.R")
## Uses eight draws of the archived full fit in archive/illustration/.
local({
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, ".here"))) {
    stop("Open GP-PU-PoE.Rproj before running this test.", call. = FALSE)
  }
  module <- file.path(root, "illustration_expert_conditioning/R")
  code <- c(
    list.files(module, "[.]R$", full.names = TRUE),
    file.path(root, "scripts/01b_expert_conditioning.R")
  )
  invisible(lapply(code, function(p) parse(file = p)))
  wf <- new.env(parent = as.environment("package:stats"))
  sys.source(file.path(module, "workflow.R"), envir = wf)
  wf$ec_dependencies()
  inputs <- wf$ec_read_inputs(root, "full")
  hashes <- inputs$fingerprints
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  n <- nrow(inputs$checkpoint$fit$draws$latent_probability)
  keep <- unique(as.integer(round(seq(1, n, length.out = 8))))
  small <- inputs
  small$checkpoint$fit$draws <- lapply(
    inputs$checkpoint$fit$draws,
    function(x) {
      if (is.matrix(x)) x[keep, , drop = FALSE] else x[keep]
    }
  )
  e <- wf$ec_environment(small)
  ## The test stops if a model fit is started.
  e$run_exact_toy_model_chain <- function(...) {
    stop("A postprocessing test tried to fit a model.")
  }
  e$run_illustration <- e$run_exact_toy_model_chain
  a <- wf$ec_compute(small, e)
  b <- wf$ec_compute(small, e)
  stopifnot(
    identical(a, b),
    nrow(a) == 200L,
    all(a$mean >= 0 & a$mean <= 1),
    all(a$lower >= 0 & a$upper <= 1 & a$lower <= a$upper),
    identical(sort(unique(a$node)), 1:100),
    identical(a$plot_id[a$expert == "cov"], match(1:100, inputs$order_nodes)),
    identical(
      a$plot_id[a$expert == "network"],
      match(1:100, inputs$order_nodes)
    )
  )
  e$expert_summary_toy <- a
  p0 <- e$plot_toy_expert(
    "cov",
    "Covariate expert",
    "dodgerblue1",
    expression(sigma(x[0]))
  )
  p1 <- e$plot_toy_expert(
    "network",
    "Network expert",
    "firebrick2",
    expression(sigma(x[1]))
  )
  stopifnot(
    identical(p0$labels$y, expression(sigma(x[0]))),
    identical(p1$labels$y, expression(sigma(x[1])))
  )
  invisible(ggplot2::ggplot_build(p0))
  invisible(ggplot2::ggplot_build(p1))
  tmp <- tempfile("expert_conditioning_test_", fileext = ".rds")
  on.exit(unlink(tmp, expand = FALSE), add = TRUE)
  saveRDS(a, tmp)
  stopifnot(
    identical(a, readRDS(tmp)),
    identical(old_kind, RNGkind()),
    identical(
      had_seed,
      exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    )
  )
  if (had_seed) {
    stopifnot(identical(
      old_seed,
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    ))
  }
  stopifnot(identical(hashes, wf$ec_fingerprints(inputs$paths)))
  ## Gaussian identity for the conditional means.
  G0 <- matrix(c(2, .3, .3, 1), 2)
  G1 <- matrix(c(1, .2, .2, 1.5), 2)
  eps <- 1e-8
  r <- c(-.7, 1.1)
  R <- chol(G0 + G1 + diag(eps, 2))
  alpha <- backsolve(R, forwardsolve(t(R), r))
  m0 <- as.vector(G0 %*% alpha)
  m1 <- as.vector(G1 %*% alpha)
  stopifnot(max(abs(m0 + m1 + eps * alpha - r)) < 1e-12)
  for (G in list(G0, G1)) {
    B <- forwardsolve(t(R), G)
    V <- G - crossprod(B)
    stopifnot(
      max(abs(V - t(V))) < 1e-12,
      min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) > -1e-12
    )
  }
  cat(
    "PASS: expert scores, conditional means, order and labels of the figure,\n",
    "RNG state and unchanged inputs.\n",
    sep = ""
  )
})
