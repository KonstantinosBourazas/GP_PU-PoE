## ESS of the saved draws, as in the original script.

core_parameter_ess_rows <- list()
node_probability_ess_rows <- list()

core_row_index <- 1L
node_row_index <- 1L

for (model_spec_now in MODEL_SPECS) {
  model_id_now <- model_spec_now$model_id
  method_now <- model_spec_now$method

  model_result_now <- get_complete_model_result(
    model_id_now
  )

  draws_now <- model_result_now$fit$draws

  ## ----------------------------------------------------------
  ## Core parameters
  ## ----------------------------------------------------------

  beta0_matrix <- matrix(
    as.numeric(draws_now$beta0),
    ncol = 1L
  )

  colnames(beta0_matrix) <- "beta0"

  core_parts <- list(
    beta0_matrix,
    as.matrix(draws_now$tau),
    as.matrix(draws_now$theta)
  )

  if (isTRUE(model_spec_now$use_pu)) {
    eta_matrix <- matrix(
      as.numeric(draws_now$eta),
      ncol = 1L
    )

    colnames(eta_matrix) <- "eta"

    core_parts <- append(
      core_parts,
      list(eta_matrix),
      after = 1L
    )
  }

  core_draw_matrix <- do.call(
    cbind,
    core_parts
  )

  if (
    is.null(colnames(core_draw_matrix)) ||
      anyDuplicated(colnames(core_draw_matrix))
  ) {
    stop(
      "The core posterior-draw columns are not uniquely named for ",
      method_now,
      "."
    )
  }

  core_ess_now <- effective_size_matrix(
    core_draw_matrix
  )

  core_parameter_ess_rows[[core_row_index]] <- data.frame(
    model_id = model_id_now,
    method = method_now,
    parameter = colnames(core_draw_matrix),
    retained_draws = nrow(core_draw_matrix),
    ESS = core_ess_now,
    ESS_fraction_of_retained = core_ess_now / nrow(core_draw_matrix),
    stringsAsFactors = FALSE
  )

  core_row_index <- core_row_index + 1L

  ## ----------------------------------------------------------
  ## ESS of the plotted latent probabilities
  ## ----------------------------------------------------------

  latent_draw_matrix <- as.matrix(
    draws_now$latent_probability
  )

  if (
    !all(
      dim(latent_draw_matrix) ==
        c(
          N_RETAINED,
          N_NODES
        )
    )
  ) {
    stop(
      "Unexpected latent-probability draw dimensions for ",
      method_now,
      "."
    )
  }

  latent_ess_now <- effective_size_matrix(
    latent_draw_matrix
  )

  node_probability_ess_rows[[node_row_index]] <- data.frame(
    model_id = model_id_now,
    method = method_now,
    node = seq_len(N_NODES),
    node_group = NODE_GROUP,
    retained_draws = nrow(latent_draw_matrix),
    ESS = latent_ess_now,
    ESS_fraction_of_retained = latent_ess_now / nrow(latent_draw_matrix),
    stringsAsFactors = FALSE
  )

  node_row_index <- node_row_index + 1L
}

core_parameter_ess <- do.call(
  rbind,
  core_parameter_ess_rows
)

node_probability_ess <- do.call(
  rbind,
  node_probability_ess_rows
)

row.names(core_parameter_ess) <- NULL
row.names(node_probability_ess) <- NULL

node_ess_summary <- do.call(
  rbind,
  lapply(
    split(
      node_probability_ess,
      node_probability_ess$model_id
    ),
    function(df_now) {
      data.frame(
        model_id = df_now$model_id[1L],
        method = df_now$method[1L],
        retained_draws = df_now$retained_draws[1L],
        min_ESS = min(
          df_now$ESS,
          na.rm = TRUE
        ),
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

row.names(node_ess_summary) <- NULL

write.csv(
  core_parameter_ess,
  file.path(
    OUT_DIR,
    "toy_seed53_ESS_core_parameters.csv"
  ),
  row.names = FALSE
)

write.csv(
  node_probability_ess,
  file.path(
    OUT_DIR,
    "toy_seed53_ESS_latent_probabilities_by_node.csv"
  ),
  row.names = FALSE
)

write.csv(
  node_ess_summary,
  file.path(
    OUT_DIR,
    "toy_seed53_ESS_latent_probabilities_summary.csv"
  ),
  row.names = FALSE
)

cat("\n================ CORE-PARAMETER ESS ================\n")

print(
  core_parameter_ess,
  row.names = FALSE
)

cat(
  "\n================ NODEWISE LATENT-PROBABILITY ESS SUMMARY ================\n"
)

print(
  node_ess_summary,
  row.names = FALSE
)

cat(
  "\nESS diagnostics were saved to:\n",
  normalizePath(
    OUT_DIR,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)
