## Four plotting functions of the original single-dataset script, unchanged.

blend_colour <- function(colour, target_colour, amount) {
  amount <- min(max(as.numeric(amount), 0), 1)
  source_rgb <- grDevices::col2rgb(colour)
  target_rgb <- grDevices::col2rgb(target_colour)
  blended_rgb <- (1 - amount) * source_rgb + amount * target_rgb
  grDevices::rgb(blended_rgb[1L, 1L], blended_rgb[2L, 1L], blended_rgb[3L, 1L], maxColorValue = 255)
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

make_text_strip <- function(label, fontsize) {
  patchwork::wrap_elements(full = grid::textGrob(
    label = label, x = 0.5, y = 0.5,
    gp = grid::gpar(fontface = "bold", fontsize = fontsize)
  ))
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
