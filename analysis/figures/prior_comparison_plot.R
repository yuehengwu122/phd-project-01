# Prior Comparison Plot
# Compares Prior 1 (Horseshoe) and Prior 2 (Nonlocal) on delta parameter
# Shows distributions on both original scale (delta) and squared scale (delta^2)

# Create prior comparison plot for delta parameter
# output_path: Optional path to save the plot (PNG format)
# Returns: Creates base R plot

plot_prior_comparison_ggplot <- function(output_path = NULL) {
  
  library(ggplot2)
  library(patchwork)
  
  # Define x range
  x_vals <- seq(0.01, 10, by = 0.01)
  
  # Calculate densities for both priors
  # Prior 1: Half-t(nu=4, s=2.7)
  prior1_orig <- 2 * dt(x_vals / 2.7, df = 4) / 2.7
  prior1_var <- (2 * dt(sqrt(x_vals) / 2.7, df = 4) / 2.7) / (2 * pmax(sqrt(x_vals), 1e-12))
  
  # Prior 2: Derived from t(nu=1.15, s=7) on delta^2
  prior2_orig <- (2 * dt((x_vals^2) / 7, df = 1.15) / 7) * 2 * x_vals
  prior2_var <- 2 * dt(x_vals / 7, df = 1.15) / 7
  
  # Create data frames
  df_orig <- data.frame(
    x = rep(x_vals, 2),
    density = c(prior1_orig, prior2_orig),
    prior = rep(c("Prior 1", "Prior 2"), each = length(x_vals))
  )
  
  df_var <- data.frame(
    x = rep(x_vals, 2),
    density = c(prior1_var, prior2_var),
    prior = rep(c("Prior 1", "Prior 2"), each = length(x_vals))
  )
  
  # Plot 1: Original scale
  p1 <- ggplot(df_orig, aes(x = x, y = density, color = prior)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("Prior 1" = "green4", "Prior 2" = "red3"),
                       labels = c(expression("Prior 1: " * delta ~ " ~ t"[4]^"+" * "(s=2.7)"),
                                 expression("Prior 2: " * delta^2 ~ " ~ t"[1.15]^"+" * "(s=7)"))) +
    labs(x = "Standard Deviations", y = "Density", color = NULL) +
    xlim(0, 10) + ylim(0, 0.8) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.key.width = unit(2, "cm")
    )
  
  # Plot 2: Variance scale  
  p2 <- ggplot(df_var, aes(x = x, y = density, color = prior)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("Prior 1" = "green4", "Prior 2" = "red3"),
                       labels = c(expression("Prior 1: " * delta ~ " ~ t"[4]^"+" * "(s=2.7)"),
                                 expression("Prior 2: " * delta^2 ~ " ~ t"[1.15]^"+" * "(s=7)"))) +
    labs(x = "Variances", y = "Density", color = NULL) +
    xlim(0, 10) + ylim(0, 0.8) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.key.width = unit(2, "cm")
    )
  
  # Combine plots
  combined_plot <- p1 + p2 + 
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Save if output path provided
  if (!is.null(output_path)) {
    ggsave(output_path, combined_plot, width = 12, height = 6, dpi = 300)
  }
  
  return(combined_plot)
}

# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# Example usage (ggplot2 version):
plot_gg <- plot_prior_comparison_ggplot("analysis/figures/output/plot-prior-delta.png")
print(plot_gg)

# Quick display without saving:
# plot_prior_comparison_ggplot()