# Simulation Study Plots
# Two experiments:
# 1. Experiment 1: Varying N with fixed lambda
# 2. Experiment 2: Varying lambda with fixed N

library(dplyr)
library(ggplot2)
library(patchwork)

# ============================================================================
# EXPERIMENT 1: Effect of Sample Size (N)
# ============================================================================

# Create plots for Experiment 1 (varying N)
# data_path: Path to the CSV file containing simulation results
# output_path: Optional path to save the plot
# Returns: Combined ggplot object
plot_experiment1_by_N <- function(
  data_path = "data/simulated/results_sim1_251025.csv",
  output_path = NULL
) {
  # Load and prepare data
  df <- read.csv(data_path, row.names = NULL)

  # Clean data types
  df <- df %>%
    mutate(
      prior = factor(
        prior,
        levels = c(1, 2),
        labels = c("Horseshoe", "Nonlocal")
      ),
      N = factor(N, levels = c(30, 100, 200)),
      delta = as.numeric(as.character(delta))
    )

  # Aggregate: median over replicates for each (prior, N, delta)
  summary_df <- df %>%
    group_by(prior, N, delta) %>%
    summarise(
      med_bf10 = median(log(bf10), na.rm = TRUE),
      q5_bf10 = quantile(log(bf10), 0.05, na.rm = TRUE),
      q95_bf10 = quantile(log(bf10), 0.95, na.rm = TRUE),
      med_sdr10 = median(log(sdr10), na.rm = TRUE),
      q5_sdr10 = quantile(log(sdr10), 0.05, na.rm = TRUE),
      q95_sdr10 = quantile(log(sdr10), 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(prior, N, delta)

  # Common aesthetics
  aes_common <- aes(
    x = delta,
    color = N,
    linetype = prior,
    group = interaction(prior, N)
  )

  # Base theme with black panel border
  theme_panel_border <- theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.width = unit(2, "cm")
  )

  # Panel 1: Bridge Sampling BF10
  p_bf <- ggplot(summary_df, aes_common) +
    # geom_errorbar(aes(ymin = q5_bf10, ymax = q95_bf10),
    #               width = 0.2, lineend = "round") +
    geom_line(aes(y = med_bf10), linewidth = 0.8) +
    geom_point(aes(y = med_bf10), size = 1.8) +
    labs(
      x = expression(delta),
      y = expression("Median " * log(BF[10]^BS)),
      color = "N",
      linetype = "Prior"
    ) +
    scale_color_manual(
      values = c("30" = "red", "100" = "darkgreen", "200" = "blue")
    ) +
    scale_linetype_manual(
      values = c("Horseshoe" = "solid", "Nonlocal" = "dashed"),
      labels = c("P1(Horseshoe)", "P2(Nonlocal)")
    ) +
    theme_minimal(base_size = 12) +
    coord_cartesian(ylim = c(-3, 30)) +
    theme_panel_border

  # Panel 2: Savage-Dickey BF10
  p_sdr <- ggplot(summary_df, aes_common) +
    # geom_errorbar(aes(ymin = q5_bf10, ymax = q95_bf10),
    #               width = 0.2, lineend = "round") +
    geom_line(aes(y = med_sdr10), linewidth = 0.8) +
    geom_point(aes(y = med_sdr10), size = 1.8) +
    labs(
      x = expression(delta),
      y = expression("Median " * log(BF[10]^SDR)),
      color = "N",
      linetype = "Prior"
    ) +
    scale_color_manual(
      values = c("30" = "red", "100" = "darkgreen", "200" = "blue")
    ) +
    scale_linetype_manual(
      values = c("Horseshoe" = "solid", "Nonlocal" = "dashed"),
      labels = c("P1(Horseshoe)", "P2(Nonlocal)")
    ) +
    theme_minimal(base_size = 12) +
    coord_cartesian(ylim = c(-3, 30)) +
    theme_panel_border

  # Combine plots
  combined_plot <- p_bf +
    p_sdr +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  # Save if output path provided
  if (!is.null(output_path)) {
    ggsave(output_path, combined_plot, width = 13, height = 6, dpi = 300)
  }

  return(combined_plot)
}

# ============================================================================
# EXPERIMENT 2: Effect of Length Scale (Lambda)
# ============================================================================

# Create plots for Experiment 2 (varying lambda)
# data_path1: Path to first CSV file containing simulation results
# data_path2: Path to second CSV file containing simulation results
# output_path: Optional path to save the plot
# Returns: Combined ggplot object
plot_experiment2_by_lambda <- function(
  data_path1 = "data/simulated/results_sim1_251025.csv",
  data_path2 = "data/simulated/results_sim2_251108.csv",
  output_path = NULL
) {
  # Load and combine data
  df1 <- read.csv(data_path1, row.names = NULL)
  df2 <- read.csv(data_path2, row.names = NULL)

  # Remove extra columns and combine
  df1 <- df1[, -1]
  df2 <- df2[, -c(1, 2)]
  df <- rbind(df1, df2)

  # Filter for N=100 and clean data types
  df <- df %>%
    filter(N == 100) %>%
    mutate(
      prior = factor(
        prior,
        levels = c(1, 2),
        labels = c("Horseshoe", "Nonlocal")
      ),
      lambda = factor(round(lambda, 2), levels = c("0.15", "0.33", "0.7")),
      delta = as.numeric(as.character(delta))
    )

  # Aggregate: median over replicates for each (prior, lambda, delta)
  summary_df <- df %>%
    group_by(prior, lambda, delta) %>%
    summarise(
      med_bf10 = median(log(bf10), na.rm = TRUE),
      q5_bf10 = quantile(log(bf10), 0.05, na.rm = TRUE),
      q95_bf10 = quantile(log(bf10), 0.95, na.rm = TRUE),
      med_sdr10 = median(log(sdr10), na.rm = TRUE),
      q5_sdr10 = quantile(log(sdr10), 0.05, na.rm = TRUE),
      q95_sdr10 = quantile(log(sdr10), 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(prior, lambda, delta)

  # Common aesthetics
  aes_common <- aes(
    x = delta,
    color = lambda,
    linetype = prior,
    group = interaction(prior, lambda)
  )

  # Base theme with black panel border
  theme_panel_border <- theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.width = unit(2, "cm")
  )

  # Panel 1: Bridge Sampling BF10
  p_bf <- ggplot(summary_df, aes_common) +
    # geom_errorbar(aes(ymin = q5_bf10, ymax = q95_bf10),
    #               width = 0.2, lineend = "round") +
    geom_line(aes(y = med_bf10), linewidth = 0.8) +
    geom_point(aes(y = med_bf10), size = 1.8) +
    labs(
      x = expression(delta),
      y = expression("Median " * log(BF[10]^BS)),
      color = expression(lambda),
      linetype = "Prior"
    ) +
    scale_color_manual(
      values = c("0.15" = "red", "0.33" = "darkgreen", "0.7" = "blue"),
      labels = c("0.15", "1/3", "0.7")
    ) +
    scale_linetype_manual(
      values = c("Horseshoe" = "solid", "Nonlocal" = "dashed"),
      labels = c("P1(Horseshoe)", "P2(Nonlocal)")
    ) +
    theme_minimal(base_size = 12) +
    coord_cartesian(ylim = c(-3, 60)) +
    theme_panel_border

  # Panel 2: Savage-Dickey BF10
  p_sdr <- ggplot(summary_df, aes_common) +
    # geom_errorbar(aes(ymin = q5_bf10, ymax = q95_bf10),
    #               width = 0.2, lineend = "round") +
    geom_line(aes(y = med_sdr10), linewidth = 0.8) +
    geom_point(aes(y = med_sdr10), size = 1.8) +
    labs(
      x = expression(delta),
      y = expression("Median " * log(BF[10]^SDR)),
      color = expression(lambda),
      linetype = "Prior"
    ) +
    scale_color_manual(
      values = c("0.15" = "red", "0.33" = "darkgreen", "0.7" = "blue"),
      labels = c("0.15", "1/3", "0.7")
    ) +
    scale_linetype_manual(
      values = c("Horseshoe" = "solid", "Nonlocal" = "dashed"),
      labels = c("P1(Horseshoe)", "P2(Nonlocal)")
    ) +
    theme_minimal(base_size = 12) +
    coord_cartesian(ylim = c(-3, 60)) +
    theme_panel_border

  # Combine plots
  combined_plot <- p_bf +
    p_sdr +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  # Save if output path provided
  if (!is.null(output_path)) {
    ggsave(output_path, combined_plot, width = 13, height = 6, dpi = 300)
  }

  return(combined_plot)
}

# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# plot_experiment1_by_N()
# plot_experiment2_by_lambda()

# Example usage for Experiment 1:
plot1 <- plot_experiment1_by_N(
  "data/simulated/results_sim1_251025.csv",
  "analysis/figures/output/plot-sim-1.png"
)
# print(plot1)

# # Example usage for Experiment 2:
plot2 <- plot_experiment2_by_lambda(
  "data/simulated/results_sim1_251025.csv",
  "data/simulated/results_sim2_251108.csv",
  "analysis/figures/output/plot-sim-2.png"
)
# print(plot2)
