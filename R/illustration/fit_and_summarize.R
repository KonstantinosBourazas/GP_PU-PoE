## Fits the four models, or reuses completed fits, and computes the summaries.

model_results <- vector("list", length(MODEL_SPECS))

names(model_results) <- vapply(
  MODEL_SPECS,
  `[[`,
  character(1),
  "model_id"
)

all_models_start <- Sys.time()

for (model_index in seq_along(MODEL_SPECS)) {
  model_spec <- MODEL_SPECS[[model_index]]
  model_results[[model_spec$model_id]] <-
    run_or_load_model(model_spec)
}

all_models_end <- Sys.time()

if (
  any(
    !vapply(
      model_results,
      function(x) is.list(x) && isTRUE(x$complete),
      logical(1)
    )
  )
) {
  stop("Not all four illustration models completed successfully.")
}

## ============================================================================
## 13. COMBINE AND SAVE FIGURE-READY OUTPUTS
## ============================================================================

node_summaries_long <- do.call(
  rbind,
  lapply(model_results, `[[`, "node_summary")
)

row.names(node_summaries_long) <- NULL

model_diagnostics <- do.call(
  rbind,
  lapply(model_results, `[[`, "diagnostics")
)

row.names(model_diagnostics) <- NULL

eta_summary <- do.call(
  rbind,
  lapply(model_results, `[[`, "eta_summary")
)

row.names(eta_summary) <- NULL

prior_calibration_completed <- do.call(
  rbind,
  lapply(model_results, `[[`, "prior_calibration")
)

row.names(prior_calibration_completed) <- NULL

model_order <- vapply(
  MODEL_SPECS,
  `[[`,
  character(1),
  "method"
)

node_summaries_long$method <- factor(
  node_summaries_long$method,
  levels = model_order
)

node_summaries_long <- node_summaries_long[
  order(node_summaries_long$method, node_summaries_long$node),
  ,
  drop = FALSE
]

node_summaries_long$method <- as.character(node_summaries_long$method)

write.csv(
  node_summaries_long,
  file.path(
    OUT_DIR,
    "toy_seed53_four_models_node_summaries_HPD80.csv"
  ),
  row.names = FALSE
)

write.csv(
  model_diagnostics,
  file.path(
    OUT_DIR,
    "toy_seed53_four_models_MCMC_diagnostics.csv"
  ),
  row.names = FALSE
)

write.csv(
  eta_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_eta_posterior_HPD80.csv"
  ),
  row.names = FALSE
)

write.csv(
  prior_calibration_completed,
  file.path(
    OUT_DIR,
    "toy_seed53_four_models_active_prior_calibration.csv"
  ),
  row.names = FALSE
)

## Summaries by group, to check the pattern of the illustration. The true
## labels are not used by the fits.
group_summary_rows <- list()
group_summary_index <- 1L

for (method_name in model_order) {
  for (group_name in c(
    "true zero",
    "hidden positive",
    "observed positive"
  )) {
    rows_now <- node_summaries_long[
      node_summaries_long$method == method_name &
        node_summaries_long$node_group == group_name,
      ,
      drop = FALSE
    ]

    group_summary_rows[[group_summary_index]] <- data.frame(
      method = method_name,
      node_group = group_name,
      n_nodes = nrow(rows_now),
      mean_posterior_predictive_mean = mean(
        rows_now$posterior_predictive_mean
      ),
      mean_HPD80_width = mean(rows_now$HPD80_width),
      median_HPD80_width = median(rows_now$HPD80_width),
      mean_conditional_T_probability = mean(
        rows_now$conditional_T_probability_mean
      ),
      mean_conditional_T_HPD80_width = mean(
        rows_now$conditional_T_HPD80_width
      ),
      stringsAsFactors = FALSE
    )

    group_summary_index <- group_summary_index + 1L
  }
}

group_summary <- do.call(rbind, group_summary_rows)
row.names(group_summary) <- NULL

write.csv(
  group_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_four_models_group_summary.csv"
  ),
  row.names = FALSE
)

## Wide table, read directly by the 2 x 2 figure.
wide_parts <- lapply(MODEL_SPECS, function(model_spec) {
  rows <- node_summaries_long[
    node_summaries_long$model_id == model_spec$model_id,
    c(
      "node",
      "posterior_predictive_mean",
      "HPD80_lower",
      "HPD80_upper",
      "HPD80_width",
      "conditional_T_probability_mean",
      "conditional_T_HPD80_lower",
      "conditional_T_HPD80_upper"
    ),
    drop = FALSE
  ]

  suffix <- model_spec$model_id
  names(rows)[-1L] <- paste0(names(rows)[-1L], "__", suffix)
  rows
})

node_summaries_wide <- data.frame(
  node = seq_len(N_NODES),
  T = T_TRUE,
  Y = Y,
  node_group = NODE_GROUP,
  stringsAsFactors = FALSE
)

for (part in wide_parts) {
  node_summaries_wide <- merge(
    node_summaries_wide,
    part,
    by = "node",
    all.x = TRUE,
    sort = FALSE
  )
}

node_summaries_wide <- node_summaries_wide[
  order(node_summaries_wide$node),
  ,
  drop = FALSE
]

write.csv(
  node_summaries_wide,
  file.path(
    OUT_DIR,
    "toy_seed53_four_models_figure_ready_wide.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 14. SAVE COMPLETE RESULT BUNDLE
## ============================================================================

total_runtime_sec <- as.numeric(difftime(
  all_models_end,
  all_models_start,
  units = "secs"
))

model_result_summaries <- lapply(
  model_results,
  function(result) {
    list(
      complete = result$complete,
      model_signature = result$model_signature,
      model_spec = result$model_spec,
      node_summary = result$node_summary,
      eta_summary = result$eta_summary,
      diagnostics = result$diagnostics,
      prior_calibration = result$prior_calibration,
      completed_at = result$completed_at,
      checkpoint_file = file.path(
        MODEL_CHECKPOINT_DIR,
        paste0(result$model_spec$model_id, "_complete.rds")
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
  n_nodes = N_NODES,
  true_zero_idx = TRUE_ZERO_IDX,
  hidden_positive_idx = HIDDEN_POS_IDX,
  observed_positive_idx = OBSERVED_POS_IDX,
  true_positive_idx = TRUE_POS_IDX,
  covariate_dgp = list(
    n_covariates = N_COVARIATES,
    background_distribution = "N(0,1)",
    positive_distribution = "N(0.5,1)"
  ),
  network_dgp = list(
    type = "two-block Bernoulli SBM",
    block_zero_idx = TRUE_ZERO_IDX,
    block_one_idx = TRUE_POS_IDX,
    p_00 = P_00,
    p_01 = P_01,
    p_11 = P_11
  ),
  modularity = list(
    seed = MODULARITY_SEED,
    q = modularity$q,
    selected_eigenvalues = modularity$selected_eigenvalues,
    eig_tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q,
    k_max = MOD_K_MAX
  ),
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
    implementation = "coda::HPDinterval(coda::as.mcmc(draws), prob=0.80)",
    equal_tail_intervals_used = FALSE
  ),
  methods = MODEL_SPECS,
  beta0_prior = beta0_prior_for_sampler,
  eta_prior = eta_prior_for_sampler,
  eta_true = length(HIDDEN_POS_IDX) / length(TRUE_POS_IDX),
  primary_illustration_quantity = paste0(
    "Posterior mean and 80% coda HPD interval of p_i=plogis(f_i)."
  ),
  additional_PU_quantity = paste0(
    "Draw-by-draw Pr(T_i=1|Y_i,parameters); for Y=0 use ",
    "eta*p_i/{1-(1-eta)*p_i}, and for Y=1 use 1."
  ),
  resume_completed_models = RESUME_COMPLETED_MODELS,
  force_fresh_run = FORCE_FRESH_RUN,
  total_runtime_sec = total_runtime_sec
)

run_config$execution_mode <- RUN_MODE
run_config$module_version <- MODULE_VERSION
run_config$configuration <- CONFIGURATION_SNAPSHOT
run_config$environment <- ENVIRONMENT_FINGERPRINT
run_config$source_files_md5 <- SOURCE_FILE_MD5
run_config$original_source_sha256 <- ORIGINAL_SOURCE_SHA256

saveRDS(
  run_config,
  file.path(
    OUT_DIR,
    "toy_seed53_four_models_run_config.rds"
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
  MODEL_SPECS,
  model_prior_bundles,
  model_result_summaries,
  node_summaries_long,
  node_summaries_wide,
  group_summary,
  model_diagnostics,
  eta_summary,
  prior_calibration_completed,
  beta0_prior_for_sampler,
  eta_prior_for_sampler,
  run_config,
  file = file.path(
    OUT_DIR,
    "TOY_SEED53_FOUR_GP_MODELS_COMPLETE_RESULTS.RData"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    OUT_DIR,
    "sessionInfo_toy_seed53_four_GP_models.txt"
  )
)

## ============================================================================
## 15. FINAL CONSOLE OUTPUT
## ============================================================================

cat("\n================ FOUR-MODEL ILLUSTRATION COMPLETE ================\n")
cat("Toy data seed           :", TOY_SEED, "\n")
cat("Models                  :", paste(model_order, collapse = ", "), "\n")
cat("Blocked warm-up/model   :", BLOCKED_BURNIN, "\n")
cat("Joint warm-up/model     :", JOINT_BURNIN, "\n")
cat("Posterior iter/model    :", POSTERIOR_ITERATIONS, "\n")
cat("Thinning                :", THIN, "\n")
cat("Retained draws/model    :", N_RETAINED, "\n")
cat("Node intervals          : 80% HPD from coda::HPDinterval\n")
cat("Total current runtime   :", format_duration(total_runtime_sec), "\n")

cat("\n================ MCMC / ILLUSTRATION DIAGNOSTICS ================\n")
print(model_diagnostics, row.names = FALSE)

cat("\n================ ETA SUMMARY ================\n")
print(eta_summary, row.names = FALSE)

cat("\nAll outputs were written to:\n")
cat(normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE), "\n")
cat("\n================ DONE ================\n")
