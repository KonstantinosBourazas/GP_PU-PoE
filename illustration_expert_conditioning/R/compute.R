## Draws of the two experts given the saved fit. Each expert is drawn from its
## marginal conditional distribution, as in the application.
draw_toy_experts <- function(fit, X_cov, Z_network, network_geometry) {
  d <- fit$draws
  P <- as.matrix(d$latent_probability)
  if (any(!is.finite(P)) || any(P < 0 | P > 1)) {
    stop("The stored probability draws are invalid.")
  }
  F <- qlogis(P)
  if (any(P == 0 | P == 1)) {
    message(
      sum(P == 0 | P == 1),
      " stored probabilities are exactly 0/1; using finite logit limits ",
      "for these entries only (numerical tail approximation)."
    )
    F[P == 0] <- qlogis(.Machine$double.xmin)
    F[P == 1] <- qlogis(1 - .Machine$double.eps / 2)
  }
  S <- nrow(F)
  n <- ncol(F)
  D2 <- sqdist(X_cov)
  XX <- tcrossprod(X_cov)
  ZZ <- tcrossprod(Z_network)
  out <- list(cov = matrix(NA_real_, S, n), network = matrix(NA_real_, S, n))

  rng_before <- .Random.seed
  on.exit(assign(".Random.seed", rng_before, envir = .GlobalEnv))
  set.seed(456)

  draw_one <- function(G, R, alpha) {
    m <- as.vector(G %*% alpha)
    B <- forwardsolve(t(R), G)
    V <- symmetrize(G - crossprod(B))
    ee <- eigen(V, symmetric = TRUE)
    if (min(ee$values) < -1e-7 * max(1, max(diag(G)))) {
      stop(
        "The conditional covariance has a negative eigenvalue beyond rounding tolerance."
      )
    }
    as.vector(m + ee$vectors %*% (sqrt(pmax(ee$values, 0)) * rnorm(n)))
  }

  for (s in seq_len(S)) {
    th <- d$theta[s, ]
    C_cov <- rbf_kernel_from_D2(
      D2,
      sigma_f = th["sigma_cov"],
      l = th["ell_cov"]
    )
    C_net <- network_kernel_from_geom(
      network_geometry,
      ell = th["ell_network"],
      sigma = th["sigma_network"]
    )
    G_cov <- C_cov + d$tau[s, "tau_cov"]^2 * XX
    G_net <- C_net + d$tau[s, "tau_network"]^2 * ZZ
    R <- chol(symmetrize(G_cov + G_net) + diag(JITTER, n))
    residual <- F[s, ] - d$beta0[s]
    alpha <- backsolve(R, forwardsolve(t(R), residual))
    out$cov[s, ] <- draw_one(G_cov, R, alpha)
    out$network[s, ] <- draw_one(G_net, R, alpha)
    if (s %% 1000L == 0L || s == S) cat("Expert conditioning:", s, "/", S, "\n")
  }
  ## Marginal draws, as in the application, not a joint sample of the two experts.
  out
}
