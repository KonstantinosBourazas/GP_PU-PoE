## Figure of the four models, from the stored summaries.

required_plot_packages <- c(
  "ggplot2",
  "patchwork"
)

missing_plot_packages <- required_plot_packages[
  !vapply(
    required_plot_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_plot_packages) > 0L) {
  stop(
    "Install the required plotting package(s): ",
    paste(missing_plot_packages, collapse = ", ")
  )
}

if (!exists("node_summaries_long")) {
  stop(
    "Object 'node_summaries_long' was not found. ",
    "Run the fit or load the saved results first."
  )
}

out_pdf_2x2_toy_HPD80 <- file.path(
  OUT_DIR,
  paste0(
    "HPD_2x2_FOUR_GP_MODELS_TOY_SEED53_",
    "SORTED_BY_PUPOE_LATENT_MEAN_80.pdf"
  )
)

## ------------------------------------------------------------
## Colours of the four models
## ------------------------------------------------------------

cols_toy_HPD80 <- c(
  "Covariates model" = "dodgerblue1",
  "Network model" = "firebrick2",
  "PoE model" = "mediumorchid2",
  "PU-PoE model" = "green3"
)

## ------------------------------------------------------------
## Construct the common node ordering
##
## Preserve the three blocks:
##   1:70   true zeros
##   71:80  hidden positives
##   81:100 observed positives
##
## Within each block, sort by the PU-PoE latent posterior mean.
## ------------------------------------------------------------

pupoe_order_reference <- node_summaries_long[
  node_summaries_long$model_id == "GP_PU_POE",
  c(
    "node",
    "node_group",
    "latent_probability_mean"
  ),
  drop = FALSE
]

if (nrow(pupoe_order_reference) != N_NODES) {
  stop(
    "Expected ",
    N_NODES,
    " GP-PU-PoE rows when constructing the plotting order."
  )
}

group_levels_toy <- c(
  "true zero",
  "hidden positive",
  "observed positive"
)

pupoe_order_reference$group_order <- match(
  pupoe_order_reference$node_group,
  group_levels_toy
)

if (anyNA(pupoe_order_reference$group_order)) {
  stop("Unexpected node_group value in the PU-PoE summaries.")
}

pupoe_order_reference <- pupoe_order_reference[
  order(
    pupoe_order_reference$group_order,
    pupoe_order_reference$latent_probability_mean,
    pupoe_order_reference$node
  ),
  ,
  drop = FALSE
]

toy_common_plot_order <- pupoe_order_reference$node

if (
  length(toy_common_plot_order) != N_NODES ||
    anyDuplicated(toy_common_plot_order) ||
    !setequal(toy_common_plot_order, seq_len(N_NODES))
) {
  stop("The common plotting order is invalid.")
}

write.csv(
  data.frame(
    plot_position = seq_len(N_NODES),
    original_node = toy_common_plot_order,
    node_group = pupoe_order_reference$node_group,
    PU_PoE_latent_posterior_mean = pupoe_order_reference$latent_probability_mean,
    stringsAsFactors = FALSE
  ),
  file.path(
    OUT_DIR,
    "toy_seed53_common_plot_order_by_PU_PoE_latent_mean.csv"
  ),
  row.names = FALSE
)

## ------------------------------------------------------------
## Block separators
## ------------------------------------------------------------

toy_block_lines_HPD80 <- data.frame(
  x = c(
    70.5,
    80.5
  )
)

## ------------------------------------------------------------
## Common plotting function
## ------------------------------------------------------------

make_hpd_plot_toy_HPD80 <- function(
  model_id_now,
  colour_now,
  title_now
) {
  df <- node_summaries_long[
    node_summaries_long$model_id == model_id_now,
    ,
    drop = FALSE
  ]

  if (nrow(df) != N_NODES) {
    stop(
      "Expected ",
      N_NODES,
      " rows for ",
      model_id_now,
      ", but found ",
      nrow(df),
      "."
    )
  }

  ## The same node order in every panel.
  df$plot_id <- match(
    df$node,
    toy_common_plot_order
  )

  df <- df[
    order(df$plot_id),
    ,
    drop = FALSE
  ]

  ## Common inferential target for all four models.
  df$plot_mean <- df$latent_probability_mean
  df$plot_lower <- df$latent_HPD80_lower
  df$plot_upper <- df$latent_HPD80_upper

  if (
    any(!is.finite(df$plot_mean)) ||
      any(!is.finite(df$plot_lower)) ||
      any(!is.finite(df$plot_upper))
  ) {
    stop(
      "Non-finite posterior summaries were found for ",
      model_id_now,
      "."
    )
  }

  if (any(df$plot_lower > df$plot_upper)) {
    stop(
      "Invalid HPD limits were found for ",
      model_id_now,
      "."
    )
  }

  mean_width_now <- mean(
    df$plot_upper -
      df$plot_lower
  )

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = plot_id
    )
  ) +

    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = plot_lower,
        ymax = plot_upper
      ),
      fill = colour_now,
      alpha = 0.20
    ) +

    ggplot2::geom_line(
      ggplot2::aes(
        y = plot_lower
      ),
      colour = colour_now,
      linewidth = 0.5
    ) +

    ggplot2::geom_line(
      ggplot2::aes(
        y = plot_upper
      ),
      colour = colour_now,
      linewidth = 0.5
    ) +

    ggplot2::geom_line(
      ggplot2::aes(
        y = plot_mean
      ),
      colour = colour_now,
      linewidth = 0.8,
      linetype = 2
    ) +

    ggplot2::geom_segment(
      data = toy_block_lines_HPD80,
      ggplot2::aes(
        x = x,
        xend = x,
        y = 0,
        yend = 1
      ),
      inherit.aes = FALSE,
      linetype = "longdash",
      linewidth = 1,
      colour = "grey35"
    ) +

    ggplot2::annotate(
      geom = "label",
      x = N_NODES / 2,
      y = 0.98,
      label = sprintf(
        "mean width = %.3f",
        mean_width_now
      ),
      vjust = 1,
      hjust = 0.5,
      size = 3.0,
      fill = "white",
      alpha = 0.90,
      label.size = 0
    ) +

    ggplot2::coord_cartesian(
      ylim = c(0, 1),
      clip = "off"
    ) +

    ggplot2::labs(
      title = title_now,
      x = "Unit, ordered within group",
      y = "Posterior probability of T = 1"
    ) +

    ggplot2::theme_minimal(
      base_size = 12
    ) +

    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 17,
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(
        size = 13
      ),
      axis.text = ggplot2::element_text(
        size = 11
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      plot.margin = ggplot2::margin(
        10,
        10,
        10,
        10
      )
    )
}

## ------------------------------------------------------------
## Four panels
## ------------------------------------------------------------

p_covariates_toy_HPD80 <- make_hpd_plot_toy_HPD80(
  model_id_now = "GP_COV",
  colour_now = cols_toy_HPD80["Covariates model"],
  title_now = "Covariates model"
)

p_network_toy_HPD80 <- make_hpd_plot_toy_HPD80(
  model_id_now = "GP_NET",
  colour_now = cols_toy_HPD80["Network model"],
  title_now = "Network model"
)

p_poe_toy_HPD80 <- make_hpd_plot_toy_HPD80(
  model_id_now = "GP_POE",
  colour_now = cols_toy_HPD80["PoE model"],
  title_now = "PoE model"
)

p_pupoe_toy_HPD80 <- make_hpd_plot_toy_HPD80(
  model_id_now = "GP_PU_POE",
  colour_now = cols_toy_HPD80["PU-PoE model"],
  title_now = "PU-PoE model"
)

layout_2x2_toy_HPD80 <- patchwork::wrap_plots(
  p_covariates_toy_HPD80,
  p_network_toy_HPD80,
  p_poe_toy_HPD80,
  p_pupoe_toy_HPD80,
  ncol = 2,
  byrow = TRUE
)

if (identical(RUN_MODE, "quick")) {
  out_pdf_2x2_toy_HPD80 <- sub("[.]pdf$", "_QUICK.pdf", out_pdf_2x2_toy_HPD80)
  layout_2x2_toy_HPD80 <- layout_2x2_toy_HPD80 +
    patchwork::plot_annotation(
      title = "Quick mode: not the results of the paper"
    )
}

write.csv(
  data.frame(
    model_id = c("GP_COV", "GP_NET", "GP_POE", "GP_PU_POE"),
    execution_mode = RUN_MODE,
    retained_draws = N_RETAINED
  ),
  file.path(OUT_DIR, "figure_run_mode.csv"),
  row.names = FALSE
)

grDevices::pdf(
  file = out_pdf_2x2_toy_HPD80,
  width = 14,
  height = 10,
  useDingbats = FALSE
)

tryCatch(print(layout_2x2_toy_HPD80), finally = grDevices::dev.off())

message(
  "Saved the figure of the four models to: ",
  out_pdf_2x2_toy_HPD80
)

## ------------------------------------------------------------
## Print eta posterior summary
## ------------------------------------------------------------

eta_report <- eta_summary[
  eta_summary$model_id == "GP_PU_POE",
  c(
    "eta_true",
    "eta_posterior_mean",
    "eta_posterior_sd",
    "eta_HPD80_lower",
    "eta_HPD80_upper"
  ),
  drop = FALSE
]

cat("\n================ ETA POSTERIOR SUMMARY ================\n")
print(eta_report, row.names = FALSE)

cat(sprintf(
  paste0(
    "\nPosterior mean eta = %.3f, posterior SD = %.3f, ",
    "80%% HPD interval = [%.3f, %.3f].\n"
  ),
  eta_report$eta_posterior_mean,
  eta_report$eta_posterior_sd,
  eta_report$eta_HPD80_lower,
  eta_report$eta_HPD80_upper
))
