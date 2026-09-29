## Summaries, ESS and plot order of the 13 fits, as in the original script.

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
        prior_role = sensitivity_results[[fit_id_now]]$fit_spec$prior_role,
        eta_prior_alpha = eta_prior_now$a_eta,
        eta_prior_beta = eta_prior_now$b_eta,
        eta_prior_mean = eta_prior_now$a_eta /
          (eta_prior_now$a_eta +
            eta_prior_now$b_eta),
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

for (row_index in seq_len(
  nrow(eta_prior_sensitivity_summary)
)) {
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
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_node_summaries_HPD80.csv"
  ),
  row.names = FALSE
)

write.csv(
  eta_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_eta_HPD80.csv"
  ),
  row.names = FALSE
)

write.csv(
  fit_diagnostics,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_MCMC_diagnostics.csv"
  ),
  row.names = FALSE
)

write.csv(
  prior_calibration_completed,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_active_priors.csv"
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

## ordered_top_nodes is loaded from source_functions.R.

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

## summarize_parameter_draws is loaded from source_functions.R.

parameter_summary <- do.call(
  rbind,
  lapply(sensitivity_results, summarize_parameter_draws)
)

row.names(parameter_summary) <- NULL

write.csv(
  sensitivity_fit_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_fit_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  topk_overlap_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_topk_overlap.csv"
  ),
  row.names = FALSE
)

write.csv(
  parameter_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_parameter_summary.csv"
  ),
  row.names = FALSE
)

## ============================================================================
## 13A. EFFECTIVE SAMPLE SIZE DIAGNOSTICS
## ============================================================================

## effective_size_matrix is loaded from source_functions.R.

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
    ESS_fraction_of_retained = core_ess_now / nrow(core_draw_matrix),
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
    ESS_fraction_of_retained = latent_ess_now / nrow(latent_draw_matrix),
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
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_ESS_core_parameters.csv"
  ),
  row.names = FALSE
)

write.csv(
  node_probability_ess,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_ESS_latent_probabilities_by_node.csv"
  ),
  row.names = FALSE
)

write.csv(
  node_ess_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_ESS_latent_probabilities_summary.csv"
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
    default_latent_posterior_mean = baseline_order_reference$latent_probability_mean,
    stringsAsFactors = FALSE
  ),
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_13fits_common_order.csv"
  ),
  row.names = FALSE
)


## Aliases used by the comparison of the 13 fits and by the 6 x 3 figure.
sensitivity_results_all <- sensitivity_results
eta_summary_all <- eta_summary
sensitivity_fit_summary_all <- sensitivity_fit_summary
node_summaries_long_all <- node_summaries_long

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
        eta_prior_mean = eta_prior_now$a_eta /
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


write.csv(
  joint_eta_comparison,
  file.path(OUT_DIR, "toy_seed53_PU_PoE_JOINT_PRIOR_STRESS_TEST_eta_HPD80.csv"),
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

cat("\n================ ALL 13 FIT-LEVEL ROBUSTNESS RESULTS ================\n")
print(sensitivity_fit_summary, row.names = FALSE)
cat("\n================ ALL 13 NODEWISE ESS SUMMARIES ================\n")
print(node_ess_summary, row.names = FALSE)
