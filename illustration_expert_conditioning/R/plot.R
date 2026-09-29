## Figure of the two experts.
plot_toy_expert <- function(expert_name, title_text, colour_now, y_label) {
  df <- expert_summary_toy[expert_summary_toy$expert == expert_name, ]
  df <- df[order(df$plot_id), ]
  ggplot2::ggplot(df, ggplot2::aes(x = plot_id)) +
    ggplot2::geom_hline(
      yintercept = 0.5,
      colour = "grey60",
      linetype = "dotted"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = colour_now,
      alpha = 0.20
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = lower),
      colour = colour_now,
      linewidth = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = upper),
      colour = colour_now,
      linewidth = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = mean),
      colour = colour_now,
      linewidth = 0.8,
      linetype = 2
    ) +
    ggplot2::geom_segment(
      data = toy_block_lines_HPD80,
      ggplot2::aes(x = x, xend = x, y = 0, yend = 1),
      inherit.aes = FALSE,
      linetype = "longdash",
      linewidth = 1,
      colour = "grey35"
    ) +
    ggplot2::annotate(
      "label",
      x = N_NODES / 2,
      y = 0.98,
      label = sprintf("mean width = %.3f", mean(df$upper - df$lower)),
      vjust = 1,
      hjust = 0.5,
      size = 3,
      fill = "white",
      alpha = 0.90,
      label.size = 0
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off") +
    ggplot2::labs(
      title = title_text,
      x = "Unit, ordered within group",
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 17, hjust = 0.5),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}


make_expert_figure <- function() {
  patchwork::wrap_plots(
    plot_toy_expert(
      "cov",
      "Covariate expert",
      cols_toy_HPD80["Covariates model"],
      expression(sigma(x[0]))
    ),
    plot_toy_expert(
      "network",
      "Network expert",
      cols_toy_HPD80["Network model"],
      expression(sigma(x[1]))
    ),
    ncol = 2
  )
}
