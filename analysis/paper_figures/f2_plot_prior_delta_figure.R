source("src/moment_prior.R")

library(ggplot2)

output_dir <- "analysis/prior_delta_figure"
output_dir <- "analysis/paper_figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

orders <- 0:3
labels <- c("d = 0", "d = 1", "d = 2", "d = 3")
colors <- c("d = 0" = "#E69F00", "d = 1" = "#56B4E9", "d = 2" = "#009E73", "d = 3" = "#CC79A7")
linetypes <- c("d = 0" = "solid", "d = 1" = "solid", "d = 2" = "solid", "d = 3" = "solid")

delta_grid <- seq(0, 12.5, length.out = 2000)

density_df <- do.call(
  rbind,
  lapply(seq_along(orders), function(i) {
    data.frame(
      delta = delta_grid,
      density = prior_moment_density(
        delta_grid,
        d = orders[i],
        kappa = get_kappa(orders[i])
      ),
      prior = factor(labels[i], levels = labels)
    )
  })
)

p <- ggplot(density_df, aes(x = delta, y = density, color = prior, linetype = prior)) +
  geom_line(linewidth = 1.05, lineend = "round") +
  scale_color_manual(values = colors, name = NULL) +
  scale_linetype_manual(values = linetypes, name = NULL) +
  scale_x_continuous(
    limits = c(0, 12.5),
    breaks = seq(0, 12, by = 2),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, 0.43),
    breaks = seq(0, 0.4, by = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = expression(delta),
    y = "Density"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.86, 0.84),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3),
    legend.key.width = unit(1.7, "lines"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

ggsave(
  filename = file.path(output_dir, "plot-prior-delta.png"),
  plot = p,
  width = 7.2,
  height = 4.6,
  dpi = 320,
  bg = "white"
)
