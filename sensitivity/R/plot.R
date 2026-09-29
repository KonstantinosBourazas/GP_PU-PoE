## Figure of the 13 fits (6 x 3 panels), with the colours, layout and order of
## the original script.
stopifnot(length(common_plot_order) == 100L, !anyDuplicated(common_plot_order),
          nrow(node_summaries_long_all) == 1300L,
          length(unique(node_summaries_long_all$fit_id)) == 13L)
N_NODES <- length(common_plot_order)
blend_colour <- function(colour, target_colour, amount) {
  amount <- min(max(as.numeric(amount), 0), 1)
  source_rgb <- grDevices::col2rgb(colour)
  target_rgb <- grDevices::col2rgb(target_colour)
  blended_rgb <- (1 - amount) * source_rgb + amount * target_rgb
  grDevices::rgb(blended_rgb[1L, 1L], blended_rgb[2L, 1L], blended_rgb[3L, 1L], maxColorValue = 255)
}

BASE_PUPOE_COLOUR <- "green3"
conservative_base <- blend_colour(BASE_PUPOE_COLOUR, "white", 0.05)
diffused_base <- blend_colour(BASE_PUPOE_COLOUR, "black", 0.05)
column_colours <- c(
  conservative = blend_colour(conservative_base, "yellow", 0.6),
  default = BASE_PUPOE_COLOUR,
  diffuse = blend_colour(diffused_base, "blue", 0.5)
)

## --------------------------------------------------------------------------
## Labels and panel map
## --------------------------------------------------------------------------

baseline_eta_beta <- sensitivity_prior_bundles$DEFAULT$eta_prior$b_eta
baseline_eta_beta_display <- as.numeric(formatC(baseline_eta_beta, format = "f", digits = 2L))
column_order <- c("conservative", "default", "diffuse")
column_headers <- c(
  conservative = "Conservative prior",
  default = "Default prior setting",
  diffuse = "Diffuse prior"
)

panel_specs <- list(
  sigma = list(
    row_label = expression(sigma[j]),
    fit_ids = c("SIGMA_CONSERVATIVE", "DEFAULT", "SIGMA_DIFFUSE"),
    titles = list("Total GP SD = 1", "Total GP SD = 2", "Total GP SD = 4")
  ),
  ell = list(
    row_label = expression(l[j]),
    fit_ids = c("ELL_CONSERVATIVE", "DEFAULT", "ELL_DIFFUSE"),
    titles = list(
      expression("Correlation" %in% "[0.05, 0.25]"),
      expression("Correlation" %in% "[0.10, 0.50]"),
      expression("Correlation" %in% "[0.10, 0.90]")
    )
  ),
  gamma = list(
    row_label = expression(gamma),
    fit_ids = c("GAMMA_CONSERVATIVE", "DEFAULT", "GAMMA_DIFFUSE"),
    titles = list(expression(0.5 * s[gamma]), expression(s[gamma]), expression(2 * s[gamma]))
  ),
  beta0 = list(
    row_label = expression(beta[0]),
    fit_ids = c("BETA0_CONSERVATIVE", "DEFAULT", "BETA0_DIFFUSE"),
    titles = list(
      expression(logit^{-1}*(beta[0]) %in% "[0.005, 0.05]"),
      expression(logit^{-1}*(beta[0]) %in% "[0.01, 0.10]"),
      expression(logit^{-1}*(beta[0]) %in% "[0.02, 0.20]")
    )
  ),
  eta = list(
    row_label = expression(eta),
    fit_ids = c("ETA_CONSERVATIVE", "DEFAULT", "ETA_DIFFUSE"),
    titles = list(
      expression(eta * " ~ " * Beta(1, 5)),
      bquote(eta * " ~ " * Beta(2, .(baseline_eta_beta_display))),
      expression(eta * " ~ " * Beta(1, 1))
    )
  )
)
row_order <- c("sigma", "ell", "gamma", "beta0", "eta")

## A compact map of the panels, saved later as in the original script.
plot_grid_map <- do.call(rbind, lapply(row_order, function(row_key_now) {
  row_spec_now <- panel_specs[[row_key_now]]
  data.frame(
    row_key = row_key_now,
    row_label = paste(deparse(row_spec_now$row_label), collapse = ""),
    column_key = column_order,
    fit_id = row_spec_now$fit_ids,
    panel_title = vapply(row_spec_now$titles, function(x) paste(deparse(x), collapse = ""), character(1)),
    stringsAsFactors = FALSE
  )
}))


make_text_strip <- function(label, fontsize) {
  patchwork::wrap_elements(full = grid::textGrob(
    label = label, x = 0.5, y = 0.5,
    gp = grid::gpar(fontface = "bold", fontsize = fontsize)
  ))
}

make_sensitivity_panel_6x3 <- function(
    fit_id_now,
    column_key_now,
    panel_title_now,
    show_x_title = FALSE,
    show_y_title = FALSE
) {
  df <- node_summaries_long_all[
    node_summaries_long_all$fit_id == fit_id_now,
    ,
    drop = FALSE
  ]
  
  if (nrow(df) != N_NODES) {
    stop(
      "Expected ",
      N_NODES,
      " rows for ",
      fit_id_now,
      "."
    )
  }
  
  df$plot_id <- match(
    df$node,
    common_plot_order
  )
  
  df <- df[
    order(df$plot_id),
    ,
    drop = FALSE
  ]
  
  df$plot_mean <- df$latent_probability_mean
  df$plot_lower <- df$latent_HPD80_lower
  df$plot_upper <- df$latent_HPD80_upper
  
  if (
    any(!is.finite(df$plot_mean)) ||
    any(!is.finite(df$plot_lower)) ||
    any(!is.finite(df$plot_upper)) ||
    any(df$plot_lower > df$plot_upper)
  ) {
    stop(
      "Invalid posterior summaries for ",
      fit_id_now,
      "."
    )
  }
  
  colour_now <- unname(
    column_colours[column_key_now]
  )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = plot_id)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = plot_lower,
        ymax = plot_upper
      ),
      fill = colour_now,
      alpha = 0.22
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = plot_lower),
      colour = colour_now,
      linewidth = 0.35
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = plot_upper),
      colour = colour_now,
      linewidth = 0.35
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = plot_mean),
      colour = colour_now,
      linewidth = 0.65,
      linetype = 2
    ) +
    ggplot2::geom_vline(
      xintercept = c(70.5, 80.5),
      linetype = "dashed",
      linewidth = 0.33,
      colour = "grey45"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      limits = c(0, 100),
      expand = ggplot2::expansion(
        mult = c(0, 0)
      )
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      limits = c(0, 1),
      expand = ggplot2::expansion(
        mult = c(0.01, 0.02)
      )
    ) +
    ggplot2::labs(
      title = panel_title_now,
      x = if (show_x_title) {
        "Unit, ordered within group"
      } else {
        NULL
      },
      y = if (show_y_title) {
        "Post. Prob. of T = 1"
      } else {
        NULL
      }
    ) +
    ggplot2::theme_minimal(
      base_size = 8.5
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "plain",
        size = 8.3,
        hjust = 0.5,
        margin = ggplot2::margin(b = 2)
      ),
      axis.title = ggplot2::element_text(
        size = 7.5
      ),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 3)
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 3)
      ),
      axis.text = ggplot2::element_text(
        size = 6.8
      ),
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.25,
        colour = "grey88"
      ),
      plot.margin = ggplot2::margin(
        4,
        4,
        4,
        4
      )
    )
}

make_vertical_text_strip <- function(
    label,
    fontsize = 11
) {
  patchwork::wrap_elements(
    full = grid::textGrob(
      label = label,
      x = 0.5,
      y = 0.5,
      rot = 90,
      gp = grid::gpar(
        fontface = "bold",
        fontsize = fontsize
      )
    )
  )
}

panel_specs_6x3 <- panel_specs

panel_specs_6x3$joint <- list(
  row_label = "Joint priors",
  fit_ids = c(
    "JOINT_CONSERVATIVE",
    "DEFAULT",
    "JOINT_DIFFUSE"
  ),
  titles = list(
    "All conservative",
    "Default",
    "All diffuse"
  )
)

row_order_6x3 <- c(
  row_order,
  "joint"
)

plot_grid_map_6x3 <- do.call(
  rbind,
  lapply(
    row_order_6x3,
    function(row_key_now) {
      row_spec_now <- panel_specs_6x3[[row_key_now]]
      
      data.frame(
        row_key = row_key_now,
        row_label = paste(
          deparse(row_spec_now$row_label),
          collapse = ""
        ),
        column_key = column_order,
        fit_id = row_spec_now$fit_ids,
        panel_title = vapply(
          row_spec_now$titles,
          function(x) paste(deparse(x), collapse = ""),
          character(1)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  plot_grid_map_6x3,
  file.path(
    OUT_DIR,
    "toy_seed53_PU_PoE_prior_sensitivity_6x3_plot_grid_map.csv"
  ),
  row.names = FALSE
)


plot_list_6x3 <- list(
  make_text_strip("", 1),
  make_text_strip(column_headers["conservative"], 12),
  make_text_strip(column_headers["default"], 12),
  make_text_strip(column_headers["diffuse"], 12)
)

for (row_index in seq_along(row_order_6x3)) {
  row_key_now <- row_order_6x3[row_index]
  row_spec_now <- panel_specs_6x3[[row_key_now]]
  
  if (identical(row_key_now, "joint")) {
    plot_list_6x3[[length(plot_list_6x3) + 1L]] <-
      make_vertical_text_strip(
        "Joint priors",
        fontsize = 11
      )
  } else {
    plot_list_6x3[[length(plot_list_6x3) + 1L]] <-
      make_text_strip(
        row_spec_now$row_label,
        15
      )
  }
  
  for (column_index in seq_along(column_order)) {
    plot_list_6x3[[length(plot_list_6x3) + 1L]] <-
      make_sensitivity_panel_6x3(
        fit_id_now =
          row_spec_now$fit_ids[column_index],
        column_key_now =
          column_order[column_index],
        panel_title_now =
          row_spec_now$titles[[column_index]],
        show_x_title =
          row_index == length(row_order_6x3),
        show_y_title =
          column_index == 1L
      )
  }
}

layout_6x3 <- patchwork::wrap_plots(
  plot_list_6x3,
  ncol = 4,
  widths = c(0.08, 1, 1, 1),
  heights = c(0.08, 1, 1, 1, 1, 1, 1),
  byrow = TRUE
) +
  patchwork::plot_annotation(
    title = if (RUN_MODE == "quick") "Quick mode: not the results of the paper" else "Prior sensitivity analysis for the PU-PoE illustration",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14,
        hjust = 0.5,
        margin = ggplot2::margin(b = 4)
      )
    )
  )

out_pdf_6x3 <- file.path(
  OUT_DIR,
  paste0("PU_PoE_PRIOR_SENSITIVITY_6x3_JOINT_STRESS_HPD80", if (RUN_MODE == "quick") "_QUICK" else "", ".pdf")
)

out_png_6x3 <- file.path(
  OUT_DIR,
  paste0("PU_PoE_PRIOR_SENSITIVITY_6x3_JOINT_STRESS_HPD80", if (RUN_MODE == "quick") "_QUICK" else "", ".png")
)

FIGURE_WIDTH_IN_6x3 <- 8.27
FIGURE_HEIGHT_IN_6x3 <- 11.69

grDevices::pdf(
  out_pdf_6x3,
  width = FIGURE_WIDTH_IN_6x3,
  height = FIGURE_HEIGHT_IN_6x3,
  useDingbats = FALSE
)

print(layout_6x3)

grDevices::dev.off()

ggplot2::ggsave(
  filename = out_png_6x3,
  plot = layout_6x3,
  width = FIGURE_WIDTH_IN_6x3,
  height = FIGURE_HEIGHT_IN_6x3,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

