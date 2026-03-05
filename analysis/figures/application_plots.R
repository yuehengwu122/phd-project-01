# =============================================================================
# Application Results Visualization
# =============================================================================
# This script creates plots for the application results, specifically posterior
# distributions of delta parameters for different fitted models.
#
# Models:
# - fit_a3_dep: Model A3 with dependent priors
# - fit_a4_dep: Model A4 with dependent priors
# =============================================================================

# Load required libraries
library(ggplot2)
library(rstan)
library(dplyr)
source("_common.R")

# =============================================================================
# Data Loading
# =============================================================================

# Load fitted models (these are result objects from fit_model function)
fit_a3_dep <- readRDS("analysis/application/fits/fit_a3_dep.rds")
fit_a4_dep <- readRDS("analysis/application/fits/fit_a4_dep.rds")

# =============================================================================
# Plotting Functions (Base R)
# =============================================================================

#' Create boxplot of posterior samples for nonlinear and linear effects
#' @param result Result object from fit_model() function
#' @param effect_type Type of effect: "nonlinear" or "linear"
#' @param sdr_values Named vector of SDR values for ranking predictors
#' @param title Plot title (optional)
#' @param ylab Y-axis label (optional)
#' @param ylim Y-axis limits as c(min, max). If NULL, auto-scales to data range.
#' @return Invisible NULL (plots directly)
plot_effect_posterior <- function(
  result,
  effect_type,
  sdr_values,
  title = NULL,
  ylab = NULL,
  ylim = NULL
) {
  fit <- result$fit
  samples <- if (effect_type == "nonlinear") {
    rstan::extract(fit)$delta
  } else if (effect_type == "linear") {
    rstan::extract(fit)$theta
  } else {
    stop("Invalid effect_type. Use 'nonlinear' or 'linear'.")
  }

  n_predictors <- ncol(samples)

  # Get predictor names from the data
  predictor_names <- colnames(result$data$X)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(n_predictors))
  }

  # Ensure SDR values are numeric
  sdr_numeric <- as.numeric(sdr_values)

  # Rank predictors by SDR values
  ranked_indices <- order(sdr_numeric, decreasing = TRUE)
  samples <- samples[, ranked_indices]
  predictor_names <- predictor_names[ranked_indices]
  sdr_ranked <- sdr_numeric[ranked_indices]

  # Set default title and ylab
  if (is.null(title)) {
    title <- if (effect_type == "nonlinear") {
      expression(paste(
        "Posterior Samples of ",
        delta[j],
        " (Nonlinear Effect)"
      ))
    } else {
      expression(paste("Posterior Samples of ", theta[j], " (Linear Effect)"))
    }
  }

  if (is.null(ylab)) {
    ylab <- if (effect_type == "nonlinear") {
      expression(delta[j])
    } else {
      expression(theta[j])
    }
  }

  # Create x-axis labels with predictor name and log(SDR) below
  x_labels <- paste0(predictor_names, "\n", sprintf("%.2f", log(sdr_ranked)))

  # Create boxplot using base R
  par(mar = c(8, 5, 4, 2))
  boxplot(
    samples,
    names = NA, # Suppress default names
    main = title,
    ylab = ylab,
    xlab = expression("Predictor with" ~ log(SDR[10])),
    col = "white",
    border = "black",
    outcol = "black",
    outpch = 3, # Use "+" for outliers
    outcex = 0.3,
    ylim = ylim,
    xaxt = "n", # Suppress x-axis
    cex.axis = 0.9
  )

  # Add custom x-axis labels without tick marks
  axis(1, at = 1:n_predictors, labels = x_labels, tick = FALSE)

  # Add horizontal line at y = 0
  abline(h = 0, lwd = 2)

  invisible(NULL)
}

# =============================================================================
# Generate Plots
# =============================================================================

# Update Generate Plots section
cat("Creating and saving individual effect posterior plots...\n")

# Compute SDR values
sdr_nonlinear_a3 <- compute_savage_dickey_ratios(fit_a3_dep)$SDR_nonlinear
sdr_linear_a3 <- compute_savage_dickey_ratios(fit_a3_dep)$SDR_linear
sdr_nonlinear_a4 <- compute_savage_dickey_ratios(fit_a4_dep)$SDR_nonlinear
sdr_linear_a4 <- compute_savage_dickey_ratios(fit_a4_dep)$SDR_linear

# Plot nonlinear effects for A3
png(
  "analysis/figures/output/delta_posterior_a3_nonlinear.png",
  width = 7,
  height = 6,
  units = "in",
  res = 300
)
plot_effect_posterior(
  fit_a3_dep,
  effect_type = "nonlinear",
  sdr_values = sdr_nonlinear_a3,
  title = expression("SDR"["10"] * "-Based Ranking of Nonlinear Effects"),
  ylab = expression("posterior of " * delta[j]),
  ylim = c(0, 10)
)
dev.off()

# Plot linear effects for A3
png(
  "analysis/figures/output/delta_posterior_a3_linear.png",
  width = 7,
  height = 6,
  units = "in",
  res = 300
)
plot_effect_posterior(
  fit_a3_dep,
  effect_type = "linear",
  sdr_values = sdr_linear_a3,
  title = expression("SDR"["10"] * "-Based Ranking of Linear Effects"),
  ylab = expression("posterior of " * theta[j]),
  ylim = c(-10, 10)
)
dev.off()

# Plot nonlinear effects for A4
png(
  "analysis/figures/output/delta_posterior_a4_nonlinear.png",
  width = 9,
  height = 7,
  units = "in",
  res = 300
)
plot_effect_posterior(
  fit_a4_dep,
  effect_type = "nonlinear",
  sdr_values = sdr_nonlinear_a4,
  title = expression("SDR"["10"] * "-Based Ranking of Nonlinear Effects"),
  ylab = expression("posterior of " * delta[j]),
  ylim = c(0, 10)
)
dev.off()

# Plot linear effects for A4
png(
  "analysis/figures/output/delta_posterior_a4_linear.png",
  width = 9,
  height = 7,
  units = "in",
  res = 300
)
plot_effect_posterior(
  fit_a4_dep,
  effect_type = "linear",
  sdr_values = sdr_linear_a4,
  title = expression("SDR"["10"] * "-Based Ranking of Linear Effects"),
  ylab = expression("posterior of " * theta[j]),
  ylim = c(-10, 10)
)
dev.off()


sdr_nonlinear_a3 <- compute_savage_dickey_ratios(
  fit_a3_dep_p2,
  delta_prior = "p2"
)$SDR_nonlinear
plot_effect_posterior(
  fit_a3_dep_p2,
  effect_type = "nonlinear",
  sdr_values = sdr_nonlinear_a3,
  title = expression("SDR"["10"] * "-Based Ranking of Nonlinear Effects"),
  ylab = expression("posterior of " * delta[j]),
  ylim = c(0, 10)
)
