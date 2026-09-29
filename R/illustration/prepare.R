## Synthetic data and prior calibration.

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

## The covariates are centered and scaled, as in the simulation and
## application code. The generated X is kept above.
X_standardization <- standardize_columns(toy$X)
X_cov <- X_standardization$scaled
colnames(X_cov) <- colnames(toy$X)
A_network <- sanitize_network_adjacency(toy$A)
Y <- as.integer(toy$Y)
T_TRUE <- as.integer(toy$T)

NODE_GROUP <- rep("true zero", N_NODES)
NODE_GROUP[HIDDEN_POS_IDX] <- "hidden positive"
NODE_GROUP[OBSERVED_POS_IDX] <- "observed positive"

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

network_geometry <- prepare_network_components_for_blocks(
  A_network,
  network_name = "toy_network"
)

write.csv(
  network_geometry$component_summary,
  file.path(OUT_DIR, "toy_seed53_network_component_summary.csv"),
  row.names = FALSE
)

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

model_prior_bundles <- lapply(
  MODEL_SPECS,
  calibrate_model_priors,
  X_cov = X_cov,
  Z_network = Z_network,
  network_geometry = network_geometry
)

names(model_prior_bundles) <- vapply(
  MODEL_SPECS,
  `[[`,
  character(1),
  "model_id"
)

prior_calibration_table <- do.call(
  rbind,
  lapply(model_prior_bundles, `[[`, "summary")
)

row.names(prior_calibration_table) <- NULL

write.csv(
  prior_calibration_table,
  file.path(OUT_DIR, "toy_seed53_prior_calibration_all_models.csv"),
  row.names = FALSE
)

cat("\n================ COMMON FIXED PRIORS ================\n")

cat(
  "beta0 ~ N(",
  beta0_prior_for_sampler$mean,
  ", ",
  beta0_prior_for_sampler$sd,
  "^2)\n",
  sep = ""
)

cat(
  "eta ~ Beta(",
  eta_prior_for_sampler$a_eta,
  ", ",
  eta_prior_for_sampler$b_eta,
  ") for GP-PU-PoE only\n",
  sep = ""
)

cat("80% node intervals     : coda::HPDinterval\n")
