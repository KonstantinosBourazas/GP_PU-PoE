## Figures of the prior sensitivity analysis over 100 datasets, with the colours
## of the original script. They average the archived summaries of each fit.

mc100_colours <- function() {
  base <- "green3"
  conservative_base <- blend_colour(base, "white", 0.05)
  diffused_base <- blend_colour(base, "black", 0.05)
  c(
    conservative = blend_colour(conservative_base, "yellow", 0.6),
    default = base,
    diffuse = blend_colour(diffused_base, "blue", 0.5)
  )
}

mc100_panels <- function() {
  list(
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
      titles = list(
        expression(0.5 * s[gamma]),
        expression(s[gamma]),
        expression(2 * s[gamma])
      )
    ),
    beta0 = list(
      row_label = expression(beta[0]),
      fit_ids = c("BETA0_CONSERVATIVE", "DEFAULT", "BETA0_DIFFUSE"),
      titles = list(
        expression(
          logit^{
            -1
          } *
            (beta[0]) %in% "[0.005, 0.05]"
        ),
        expression(
          logit^{
            -1
          } *
            (beta[0]) %in% "[0.01, 0.10]"
        ),
        expression(
          logit^{
            -1
          } *
            (beta[0]) %in% "[0.02, 0.20]"
        )
      )
    ),
    eta = list(
      row_label = expression(eta),
      fit_ids = c("ETA_CONSERVATIVE", "DEFAULT", "ETA_DIFFUSE"),
      titles = list(
        expression(eta * " ~ " * Beta(1, 5)),
        expression(eta * " ~ " * Beta(2, 6.39)),
        expression(eta * " ~ " * Beta(1, 1))
      )
    ),
    joint = list(
      row_label = "Joint priors",
      fit_ids = c("JOINT_CONSERVATIVE", "DEFAULT", "JOINT_DIFFUSE"),
      titles = list("All conservative", "Default", "All diffuse")
    )
  )
}

mc100_plot_input_check <- function(bundle) {
  mc100_assert(
    identical(bundle$format_version, "mc100_original_node_order") &&
      identical(
        bundle$prediction_index,
        "original_node_no_probability_ordering"
      ),
    "Old or incompatible combined bundle. Run ACTION 'all' once."
  )
  d <- bundle$plot_means
  et <- bundle$eta_means
  mc100_assert(
    nrow(d) == 1300L &&
      setequal(unique(d$fit_id), mc100_ids()) &&
      all(d$n_replications == 100L),
    "Expected prediction means for 13 settings and 100 positions."
  )
  for (id in mc100_ids()) {
    dd <- d[d$fit_id == id, ]
    mc100_assert(
      !anyDuplicated(dd$position) && setequal(dd$position, 1:100),
      "Bad prediction positions."
    )
  }
  vals <- unlist(d[c(
    "Mean_posterior_mean",
    "Mean_HPD80_lower",
    "Mean_HPD80_upper"
  )])
  mc100_assert(
    all(is.finite(vals)) &&
      all(vals >= 0 & vals <= 1) &&
      all(d$Mean_HPD80_lower <= d$Mean_HPD80_upper),
    "Invalid mean prediction summaries."
  )
  mc100_assert(
    nrow(et) == 13L &&
      !anyDuplicated(et$Setting) &&
      setequal(et$Setting, mc100_ids()) &&
      all(et$Replications == 100L) &&
      all(et$True_eta == 1 / 3),
    "Invalid aggregate eta table."
  )
  invisible(TRUE)
}

mc100_prediction_plot <- function(bundle) {
  d <- bundle$plot_means
  ## Units keep their indices 1:100 in every dataset and setting, without sorting.
  ## An index is a position in the simulation design, not the same unit across datasets.
  node_summaries_long_all <- data.frame(
    fit_id = d$fit_id,
    node = d$position,
    latent_probability_mean = d$Mean_posterior_mean,
    latent_HPD80_lower = d$Mean_HPD80_lower,
    latent_HPD80_upper = d$Mean_HPD80_upper
  )
  common_plot_order <- 1:100
  N_NODES <- 100L
  column_colours <- mc100_colours()
  panel <- make_sensitivity_panel_6x3
  environment(panel) <- environment() # the panel function reads the objects above
  panels <- mc100_panels()
  rows <- names(panels)
  roles <- c("conservative", "default", "diffuse")
  headers <- c("Conservative prior", "Default prior setting", "Diffuse prior")
  plots <- c(
    list(make_text_strip("", 1)),
    lapply(headers, make_text_strip, fontsize = 12)
  )
  for (j in seq_along(rows)) {
    spec <- panels[[rows[j]]]
    plots[[length(plots) + 1L]] <- if (rows[j] == "joint") {
      make_vertical_text_strip("Joint priors", 11)
    } else {
      make_text_strip(spec$row_label, 15)
    }
    for (k in seq_along(roles)) {
      p <- panel(
        fit_id_now = spec$fit_ids[k],
        column_key_now = roles[k],
        panel_title_now = spec$titles[[k]],
        show_x_title = j == length(rows),
        show_y_title = k == 1L
      )
      ## Axis labels for Monte Carlo means.
      p <- p +
        ggplot2::labs(
          x = if (j == length(rows)) "Unit" else NULL,
          y = if (k == 1L) "Mean post. prob. of T = 1" else NULL
        )
      plots[[length(plots) + 1L]] <- p
    }
  }
  patchwork::wrap_plots(
    plots,
    ncol = 4,
    widths = c(0.08, 1, 1, 1),
    heights = c(0.08, rep(1, 6)),
    byrow = TRUE
  ) +
    patchwork::plot_annotation(
      title = "Prior sensitivity analysis for the PU-PoE illustration",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 14,
          hjust = 0.5,
          margin = ggplot2::margin(b = 4)
        ),
        plot.margin = ggplot2::margin(6, 6, 6, 6)
      )
    )
}

mc100_eta_plot_data <- function(bundle) {
  panels <- mc100_panels()
  roles <- c("conservative", "default", "diffuse")
  label_keys <- c(
    sigma = 'atop("GP amplitude", sigma[j])',
    ell = 'atop("Lengthscale", l[j])',
    gamma = 'atop("Structured mean", gamma)',
    beta0 = 'atop("Intercept", beta[0])',
    eta = 'atop("False negative rate", eta)',
    joint = 'atop("Joint", "priors")'
  )
  out <- mc100_bind(lapply(names(panels), function(block) {
    ids <- panels[[block]]$fit_ids
    d <- bundle$eta_means[match(ids, bundle$eta_means$Setting), , drop = FALSE]
    data.frame(
      block = block,
      facet_label = unname(label_keys[block]),
      role = roles,
      x = 1:3,
      fit_id = ids,
      Replications = d$Replications,
      Mean_posterior_mean = d$Mean_posterior_mean,
      Mean_HPD80_lower = d$Mean_HPD80_lower,
      Mean_HPD80_upper = d$Mean_HPD80_upper,
      Coverage_HPD80 = d$Coverage_HPD80,
      stringsAsFactors = FALSE
    )
  }))
  out$role <- factor(out$role, levels = roles)
  out$facet_label <- factor(out$facet_label, levels = unname(label_keys))
  mc100_assert(
    nrow(out) == 18L &&
      length(unique(out$fit_id)) == 13L &&
      sum(out$fit_id == "DEFAULT") == 6L,
    "Eta plot must have 18 bars representing 13 distinct fits."
  )
  out
}

mc100_eta_plot <- function(bundle) {
  d <- mc100_eta_plot_data(bundle)
  colours <- mc100_colours()
  ymax <- min(1, ceiling((max(d$Mean_HPD80_upper) + 0.02) * 10) / 10)
  ymax <- max(ymax, 0.4)
  ggplot2::ggplot(d, ggplot2::aes(x = x, colour = role)) +
    ggplot2::geom_hline(
      yintercept = 1 / 3,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Mean_HPD80_lower, ymax = Mean_HPD80_upper),
      width = 0.30,
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = Mean_posterior_mean),
      size = 2.7,
      shape = 16
    ) +
    ggplot2::facet_grid(. ~ facet_label, labeller = ggplot2::label_parsed) +
    ggplot2::scale_colour_manual(
      values = colours,
      breaks = c("conservative", "default", "diffuse"),
      labels = c("Conservative prior", "Default prior", "Diffuse prior"),
      name = NULL,
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      limits = c(0.45, 3.55),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, ymax, by = 0.1),
      limits = c(0, ymax),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = expression(
        "Average posterior means and 80% HPD limits for " * eta
      ),
      x = NULL,
      y = expression(eta)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 3.2))
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      strip.text.x = ggplot2::element_text(
        size = 10.5,
        margin = ggplot2::margin(t = 5, b = 7)
      ),
      strip.background = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey88",
        linewidth = 0.25
      ),
      panel.spacing.x = grid::unit(0.65, "lines"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9.4),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(0.1, "lines"),
      axis.title.y = ggplot2::element_text(
        size = 10,
        margin = ggplot2::margin(r = 7)
      ),
      axis.text.y = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(8, 10, 7, 8)
    )
}

mc100_caption_text <- function() {
  c(
    "PREDICTION FIGURE",
    paste(
      "Prior sensitivity over 100 independently generated datasets. For each prior",
      "setting, results are averaged at the original unit indices, without sorting by",
      "posterior probabilities. Dashed curves show average posterior means; shaded",
      "bands are bounded by the averages of the lower and upper within-fit 80% HPD limits.",
      "Indices 1-70 are true zeros, 71-80 hidden positives, and 81-100 observed positives.",
      "An index is a position in the simulation design, not the same unit across datasets."
    ),
    "",
    "ETA FIGURE",
    paste(
      "Prior sensitivity of the false negative rate eta = Pr(Y=0 | T=1). Points show",
      "the average of the 100 posterior means. Vertical bars join the averages of the",
      "lower and upper within-fit 80% HPD limits. The horizontal dashed line is the true",
      "eta = 1/3. Colours distinguish conservative, default and diffuse specifications.",
      "The default result is repeated in each group; 18 bars represent 13 distinct settings."
    ),
    "",
    "INTERPRETATION",
    paste(
      "The averaged HPD limits are averages of the limits within fits. They are not a",
      "pooled posterior interval or a confidence interval. Coverage is computed from the",
      "100 intervals within fits."
    ),
    "",
    "TABLES",
    paste(
      "AUC summaries and Pearson/Spearman correlations use paired datasets.",
      "Differences are setting minus DEFAULT within each replication.",
      "MCSE_mean_difference = SD of paired differences / sqrt(100).",
      "Q05 and Q95 describe the spread over the 100 datasets.",
      "The ranges table reports both the range of the 13 setting means and the within-dataset",
      "range across the 13 settings. Node-ranking correlations compare posterior-mean scores",
      "within each dataset, separately over all units and over the 80 unlabeled units."
    ),
    "",
    "The HPD limits and ESS of each fit are the archived values."
  )
}

mc100_save_figure <- function(plot, stem, width, height) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(stem, ".pdf"),
    plot = plot,
    device = "pdf",
    useDingbats = FALSE,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    paste0(stem, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  cat("Saved: ", paste0(stem, ".pdf"), "\n", sep = "")
}

mc100_make_figures <- function(bundle, out) {
  for (p in c("ggplot2", "patchwork")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Missing plotting package: ", p)
    }
  }
  mc100_plot_input_check(bundle)
  mc100_save_figure(
    mc100_prediction_plot(bundle),
    file.path(out, "figures", "PU_PoE_PRIOR_SENSITIVITY_MC100_predictions_6x3"),
    8.27,
    11.69
  )
  mc100_save_figure(
    mc100_eta_plot(bundle),
    file.path(out, "figures", "PU_PoE_PRIOR_SENSITIVITY_MC100_eta_wide"),
    10.5,
    3.9
  )
  mc100_table(
    mc100_eta_plot_data(bundle),
    file.path(out, "combined", "eta_plot_18_bars_13_settings.csv")
  )
  colour_table <- data.frame(
    role = names(mc100_colours()),
    colour = unname(mc100_colours()),
    stringsAsFactors = FALSE
  )
  mc100_table(colour_table, file.path(out, "verification", "plot_colours.csv"))
  writeLines(
    mc100_caption_text(),
    file.path(out, "figures", "Figure_captions.txt")
  )
  invisible(TRUE)
}
