# =============================================================================
# GP Hyperparameter Effects Visualization
# =============================================================================
# This script creates visualizations showing how GP trends change with
# different values of delta (amplitude) and lambda (length scale) parameters.
#
# Grid: 4 rows (delta) × 3 columns (lambda) = 12 subplots
# Delta values: 0, 1, 2, 3
# Lambda values: 0.1, 0.3, 0.7
# =============================================================================

# Load required libraries
library(ggplot2)
library(gridExtra)
library(scales)
library(dplyr)
source("_common.R")

# =============================================================================
# GP Functions
# =============================================================================

#' Squared Exponential (SE) kernel
#' @param x1 First set of points
#' @param x2 Second set of points  
#' @param delta Amplitude parameter
#' @param lambda Length scale parameter
se_kernel <- function(x1, x2, delta, lambda) {
  # Distance matrix
  D <- outer(x1, x2, "-")
  D2 <- D^2
  
  # SE kernel: delta^2 * exp(-0.5 * D2 / lambda^2)
  delta^2 * exp(-0.5 * D2 / lambda^2)
}

#' Generate GP samples
#' @param x Input points
#' @param delta Amplitude parameter
#' @param lambda Length scale parameter
#' @param n_samples Number of samples to draw
#' @param seed Random seed
generate_gp_samples <- function(x, delta, lambda, n_samples = 5, seed = 123) {
  set.seed(seed)
  
  # Compute covariance matrix
  K <- se_kernel(x, x, delta, lambda)
  
  # Add small jitter for numerical stability
  K <- K + diag(1e-6, length(x))
  
  # Sample from multivariate normal
  L <- chol(K)
  samples <- matrix(rnorm(length(x) * n_samples), ncol = n_samples)
  samples <- t(L) %*% samples
  
  return(samples)
}

# =============================================================================
# Data Generation
# =============================================================================

# Define parameter grid
delta_values <- c(0, 1, 2, 3)
lambda_values <- c(0.1, 0.3, 0.7)

# Input points for GP
x <- seq(-0.5, 0.5, length.out = 100)

# Generate all combinations
results_list <- list()

for (i in seq_along(delta_values)) {
  for (j in seq_along(lambda_values)) {
    delta <- delta_values[i]
    lambda <- lambda_values[j]
    
    # Generate GP samples
    if (delta == 0) {
      # Special case: delta = 0 means no GP effect
      samples <- matrix(0, nrow = length(x), ncol = 5)
    } else {
      samples <- generate_gp_samples(x, delta, lambda, n_samples = 5)
    }
    
    # Create data frame
    df <- data.frame(
      x = rep(x, 5),
      y = as.vector(samples),
      sample = rep(1:5, each = length(x)),
      delta = paste0("δ = ", delta),
      lambda = paste0("λ = ", lambda),
      delta_val = delta,
      lambda_val = lambda
    )
    
    results_list[[paste0("d", i, "_l", j)]] <- df
  }
}

# Combine all results
plot_data <- do.call(rbind, results_list)

# Set factor levels for proper ordering
plot_data$delta <- factor(plot_data$delta, 
                         levels = paste0("δ = ", delta_values))
plot_data$lambda <- factor(plot_data$lambda, 
                          levels = paste0("λ = ", lambda_values))

# =============================================================================
# Visualization
# =============================================================================

#' Create GP hyperparameter effects plot
create_hyperparameter_plot <- function(data) {
  # Calculate global axis limits for consistency
  xlim_range <- range(data$x)
  ylim_range <- range(data$y)
  
  p <- ggplot(data, aes(x = x, y = y, group = sample)) +
    geom_line(color = "black", alpha = 0.7, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
    facet_grid(delta ~ lambda) +
    coord_cartesian(xlim = xlim_range, ylim = ylim_range) +
    labs(x = NULL, y = NULL) +  # Remove all labels for paper
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      strip.text.y = element_text(angle = 0),  # Keep row labels horizontal and on left
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # Add black frame
      panel.grid = element_blank(),  # Remove all grid lines
      legend.position = "none",
      panel.spacing = unit(0.8, "lines"),
      axis.text = element_text(size = 10)
    )
  
  return(p)
}

# Create the main plot
gp_effects_plot <- create_hyperparameter_plot(plot_data)

# Display the plot
print(gp_effects_plot)

# =============================================================================
# Save Plots
# =============================================================================

# Create figures directory if it doesn't exist
if (!dir.exists("analysis/figures/output")) {
  dir.create("analysis/figures/output", recursive = TRUE)
}

# Save main hyperparameter grid plot
ggsave("analysis/figures/output/gp_hyperparameter_grid.png", 
       gp_effects_plot, 
       width = 12, height = 10, dpi = 300)

# =============================================================================
# Additional Analysis: Parameter Interpretation
# =============================================================================

#' Create parameter interpretation plots
create_interpretation_plots <- function() {
  
  # Plot 1: Effect of delta (amplitude) with fixed lambda
  lambda_fixed <- 0.3
  x_demo <- seq(-0.5, 0.5, length.out = 50)
  
  delta_demo_data <- data.frame()
  for (delta in c(0.5, 1, 2, 3)) {
    samples <- generate_gp_samples(x_demo, delta, lambda_fixed, n_samples = 1, seed = 42)
    df <- data.frame(
      x = x_demo,
      y = as.vector(samples),
      delta = paste0("δ = ", delta)
    )
    delta_demo_data <- rbind(delta_demo_data, df)
  }
  
  p1 <- ggplot(delta_demo_data, aes(x = x, y = y, color = delta)) +
    geom_line(linewidth = 1.2) +
    scale_color_viridis_d(option = "viridis") +
    labs(title = "Effect of Amplitude (δ)", 
         subtitle = paste("Fixed λ =", lambda_fixed),
         x = "Input (x)", y = "f(x)") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # Plot 2: Effect of lambda (length scale) with fixed delta
  delta_fixed <- 2
  
  lambda_demo_data <- data.frame()
  for (lambda in c(0.1, 0.3, 0.7, 1.5)) {
    samples <- generate_gp_samples(x_demo, delta_fixed, lambda, n_samples = 1, seed = 42)
    df <- data.frame(
      x = x_demo,
      y = as.vector(samples),
      lambda = paste0("λ = ", lambda)
    )
    lambda_demo_data <- rbind(lambda_demo_data, df)
  }
  
  p2 <- ggplot(lambda_demo_data, aes(x = x, y = y, color = lambda)) +
    geom_line(linewidth = 1.2) +
    scale_color_viridis_d(option = "plasma") +
    labs(title = "Effect of Length Scale (λ)", 
         subtitle = paste("Fixed δ =", delta_fixed),
         x = "Input (x)", y = "f(x)") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  return(list(amplitude_effect = p1, lengthscale_effect = p2))
}

# Create interpretation plots
interpretation_plots <- create_interpretation_plots()

# Display interpretation plots
print(interpretation_plots$amplitude_effect)
print(interpretation_plots$lengthscale_effect)

