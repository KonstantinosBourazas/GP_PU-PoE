## Extracted from GP_four_model_illustrtation.R. See docs/SOURCE_MAP.csv.
## Original mathematical operations and their order are retained.

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

  ## Full mode keeps the schedule of the paper; quick mode uses a shorter one
  ## with the same transition kernels.
  validate_illustration_schedule(
    n_blocked_burnin, n_joint_burnin, n_sample, thin,
    paper_schedule = identical(RUN_MODE, "full")
  )

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
  ## Phase 1: correlated blocked warm-up (length determined by configuration)
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
  ## Phase 2: joint warm-up/tuning (length determined by configuration)
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
  ## Phase 3: frozen posterior sampling (length and thinning from configuration)
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

