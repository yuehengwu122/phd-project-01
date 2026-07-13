source("src/_setup.R")

library(ggplot2)

output_dir <- "analysis/paper_figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

style_posterior_plot <- function(plot_obj, parameter, prior_color, ylim = NULL) {
  stopifnot(inherits(plot_obj, "ggplot"))

  plot_obj$layers[[1]]$aes_params$fill <- "#6B7280"
  plot_obj$layers[[1]]$aes_params$alpha <- 0.82
  plot_obj$layers[[1]]$aes_params$color <- "white"

  if (length(plot_obj$layers) >= 2L) {
    plot_obj$layers[[2]]$aes_params$color <- prior_color
    plot_obj$layers[[2]]$aes_params$colour <- prior_color
    plot_obj$layers[[2]]$aes_params$linewidth <- 0.95
  }

  if (identical(parameter, "theta") && length(plot_obj$layers) >= 3L) {
    plot_obj$layers[[3]]$aes_params$color <- "grey45"
    plot_obj$layers[[3]]$aes_params$colour <- "grey45"
    plot_obj$layers[[3]]$aes_params$linewidth <- 0.35
    plot_obj$layers[[3]]$aes_params$linetype <- "22"
  }

  plot_obj +
    labs(
      title = NULL,
      x = if (identical(parameter, "delta")) expression(delta) else expression(theta),
      y = "Density"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = ylim) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "grey96", color = "grey70", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9.5, color = "black"),
      legend.position = "none",
      panel.spacing = grid::unit(0.9, "lines"),
      plot.margin = margin(8, 8, 6, 8)
    )
}

save_posterior_plot <- function(fit_path, plot_kind, out_name, ylim) {
  fit_full <- readRDS(fit_path)
  n_panels <- length(fit_full$predictor_names)
  n_rows <- ceiling(n_panels / 3)

  if (identical(plot_kind, "delta")) {
    plot_obj <- plot_posterior_delta(
      fit_full,
      bins = 32,
      plot_prior = TRUE,
      title = NULL
    )
    plot_obj <- style_posterior_plot(plot_obj, "delta", prior_color = "red", ylim = ylim)
  } else if (identical(plot_kind, "theta")) {
    plot_obj <- plot_posterior_theta(
      fit_full,
      bins = 32,
      plot_prior = TRUE,
      title = NULL
    )
    plot_obj <- style_posterior_plot(plot_obj, "theta", prior_color = "red", ylim = ylim)
  } else {
    stop("Unknown plot_kind: ", plot_kind)
  }

  ggsave(
    filename = file.path(output_dir, out_name),
    plot = plot_obj,
    width = 9.2,
    height = 2.45 * n_rows + 0.55,
    dpi = 320,
    bg = "white"
  )
}

save_posterior_plot(
  fit_path = "results/fit_03_diabetes_full_select.rds",
  plot_kind = "delta",
  out_name = "plot-a3-post-delta.png",
  ylim = c(0, 0.5)
)

save_posterior_plot(
  fit_path = "results/fit_03_diabetes_full_select.rds",
  plot_kind = "theta",
  out_name = "plot-a3-post-theta.png",
  ylim = c(0, 1.5)
)

save_posterior_plot(
  fit_path = "results/fit_04_boston_full_select.rds",
  plot_kind = "delta",
  out_name = "plot-a4-post-delta.png",
  ylim = c(0, 0.8)
)

save_posterior_plot(
  fit_path = "results/fit_04_boston_full_select.rds",
  plot_kind = "theta",
  out_name = "plot-a4-post-theta.png",
  ylim = c(0, 1.6)
)
