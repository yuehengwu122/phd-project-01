# Savage-Dickey Ratio Illustration
# Demonstrates the Savage-Dickey ratio method for Bayes factor calculation
# Shows prior and posterior densities with point evaluations at null hypothesis

# Create Savage-Dickey ratio illustration plot
# output_path: Optional path to save the plot (PNG format)
# Returns: Creates base R plot demonstrating SDR method
plot_savage_dickey_illustration <- function(output_path = NULL) {
  set.seed(42)

  # Create prior density: Half-Normal(0, 1)
  xgrid <- seq(0, 4, length.out = 500)
  prior_density <- 2 * dnorm(xgrid, mean = 0, sd = 1)

  # Simulate posterior draws: centered at 0.5, more concentrated
  posterior_draws <- abs(rnorm(1e4, mean = 0.65, sd = 0.3))

  # Estimate posterior density using logspline with nonnegative support
  library(logspline)
  fit_post <- logspline(posterior_draws, lbound = 0)
  posterior_density <- dlogspline(xgrid, fit_post)

  # Evaluate densities at 0
  prior_0 <- 2 * dnorm(0, 0, 1)
  post_0 <- dlogspline(0, fit_post)

  # Set up plotting parameters
  par(cex.lab = 1.2, cex.main = 1.6, cex.axis = 1, mfrow = c(1, 1))

  # Create histogram and overlay densities
  hist(
    posterior_draws,
    freq = FALSE,
    breaks = 50,
    xlim = c(-0.5, 3),
    ylim = c(0, 1.5),
    las = 1,
    border = FALSE,
    col = adjustcolor("grey", alpha.f = 0.4),
    xlab = expression(theta),
    ylab = "Density",
    main = "",
    xaxt = "n"
  )

  # Add custom x-axis
  axis(1, at = seq(-0.5, 3, by = 0.5))

  # Add density curves
  lines(xgrid, posterior_density, lwd = 2)
  lines(xgrid, prior_density, lwd = 1, col = 'blue')

  # Add points at theta = 0
  points(0, post_0, pch = 19, cex = 1.5)
  points(0, prior_0, pch = 19, cex = 1.5, col = 'blue')

  # Add vertical segments showing density evaluations
  segments(x0 = 0, y0 = 0, x1 = 0, y1 = post_0, col = "black", lty = 2, lwd = 2)
  segments(x0 = 0, y0 = 0, x1 = 0, y1 = prior_0, col = "blue", lty = 2)

  # Add labels
  text(1.3, 1.05, "Posterior", cex = 1.5)
  text(1.8, 0.25, "Prior", cex = 1.5, col = 'blue')

  text(-0.32, post_0, "post(0)", cex = 1.5)
  text(-0.32, prior_0, "prior(0)", cex = 1.5, col = 'blue')

  # Save plot if output path provided
  if (!is.null(output_path)) {
    # Save current plot to file
    dev.copy(png, filename = output_path, width = 800, height = 600, res = 150)
    dev.off()
    cat("Plot saved to:", output_path, "\n")
  }
}

# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# Example usage:
# plot_savage_dickey_illustration("analysis/figures/output/savage-dickey-illustration.png")

# Quick display without saving:
plot_savage_dickey_illustration()

# ============================================================================
# 2d savage-dickey ratio (abandoned)
# ============================================================================

# ### prior

# library(MASS)
# library(plotly)

# # Simulate prior draws
# set.seed(1)
# delta <- rnorm(1e4, 0, 1)
# xi <- rnorm(1e4, 0, 1)

# # Estimate 2D density over grid
# grid_n <- 100
# kde <- kde2d(delta, xi, n = grid_n, lims = c(0, 4, 0, 4))

# # Build surface plot
# p <- plot_ly(showscale = FALSE) %>%
#   add_surface(x = kde$x, y = kde$y, z = kde$z,
#               colorscale = list(c(0, 1), c("lightgrey", "grey")),
#               opacity = 0.9,
#               name = "Prior Density")

# # Highlight slice at delta = 0 (i.e., along xi-axis)
# delta0_idx <- which.min(abs(kde$x - 0))
# z_delta0 <- kde$z[delta0_idx, ]
# p <- p %>% add_trace(
#   x = kde$y,
#   y = rep(0, grid_n),
#   z = z_delta0,
#   type = "scatter3d",
#   mode = "lines",
#   line = list(color = "red", width = 6),
#   name = "delta = 0"
# )

# # Highlight slice at xi = 0 (i.e., along delta-axis)
# xi0_idx <- which.min(abs(kde$y - 0))
# z_xi0 <- kde$z[, xi0_idx]
# p <- p %>% add_trace(
#   x = rep(0, grid_n),
#   y = kde$x,
#   z = z_xi0,
#   type = "scatter3d",
#   mode = "lines",
#   line = list(color = "blue", width = 6),
#   name = "xi = 0"
# )

# # Add vertical lines at fixed xi along delta = 0
# xi_sample_idx <- seq(1, grid_n, length.out = 50)
# for (i in xi_sample_idx) {
#   p <- p %>% add_trace(
#     x = c(kde$y[i], kde$y[i]),
#     y = c(0, 0),
#     z = c(0, kde$z[delta0_idx, i]),
#     type = "scatter3d",
#     mode = "lines",
#     line = list(color = "red", width = 2),
#     showlegend = FALSE
#   )
# }

# # Add vertical lines at fixed delta along xi = 0
# delta_sample_idx <- seq(1, grid_n, length.out = 50)
# for (j in delta_sample_idx) {
#   p <- p %>% add_trace(
#     x = c(0, 0),
#     y = c(kde$x[j], kde$x[j]),
#     z = c(0, kde$z[j, xi0_idx]),
#     type = "scatter3d",
#     mode = "lines",
#     line = list(color = "blue", width = 2),
#     showlegend = FALSE
#   )
# }

# # Final layout
# p <- p %>% layout(
#   scene = list(
#     xaxis = list(title = "delta", showgrid = FALSE),
#     yaxis = list(title = "xi", showgrid = FALSE),
#     zaxis = list(title = "Density", showgrid = FALSE, range = c(0, 1)),
#     aspectratio = list(x = 1, y = 1, z = 0.7)
#   ),
#   title = "Joint Prior Density"
# )

# ### posterior

# # Simulate prior draws
# delta <-abs(rnorm(1e4, 1, 0.4))
# xi <- abs(rnorm(1e4, 0, 0.03))

# # delta <- fitGP_post$xi_1[,1, 1]
# # xi <- fitGP_post$xi_1[,1, 4]

# # Estimate 2D density over grid
# grid_n <- 100
# kde <- kde2d(delta, xi, n = grid_n, lims = c(0, 6, 0, 6))

# # Build surface plot
# p <- plot_ly(showscale = FALSE) %>%
#   add_surface(x = kde$x, y = kde$y, z = t(kde$z),
#               colorscale = list(c(0, 1), c("lightgrey", "grey")),
#               opacity = 0.9,
#               name = "Prior Density")

# # Highlight slice at delta = 0 (i.e., along xi-axis)
# delta0_idx <- which.min(abs(kde$x - 0))
# z_delta0 <- kde$z[delta0_idx, ]
# p <- p %>% add_trace(
#   x = rep(0, grid_n),
#   y = kde$y,
#   z = z_delta0,
#   type = "scatter3d",
#   mode = "lines",
#   line = list(color = "red", width = 6),
#   name = "delta = 0"
# )

# # Highlight slice at xi = 0 (i.e., along delta-axis)
# xi0_idx <- which.min(abs(kde$y - 0))
# z_xi0 <- kde$z[, xi0_idx]
# p <- p %>% add_trace(
#   x = kde$x,
#   y = rep(0, grid_n),
#   z = z_xi0,
#   type = "scatter3d",
#   mode = "lines",
#   line = list(color = "blue", width = 6),
#   name = "xi = 0"
# )

# # Add vertical lines at fixed xi along delta = 0
# xi_sample_idx <- seq(1, grid_n, length.out = 50)
# for (i in xi_sample_idx) {
#   p <- p %>% add_trace(
#     x = c(0, 0),
#     y = c(kde$y[i], kde$y[i]),
#     z = c(0, kde$z[delta0_idx, i]),
#     type = "scatter3d",
#     mode = "lines",
#     line = list(color = "red", width = 2),
#     showlegend = FALSE
#   )
# }

# # Add vertical lines at fixed delta along xi = 0
# delta_sample_idx <- seq(1, grid_n, length.out = 50)
# for (j in delta_sample_idx) {
#   p <- p %>% add_trace(
#     x = c(kde$x[j], kde$x[j]),
#     y = c(0, 0),
#     z = c(0, kde$z[j, xi0_idx]),
#     type = "scatter3d",
#     mode = "lines",
#     line = list(color = "blue", width = 2),
#     showlegend = FALSE
#   )
# }

# # Final layout
# p <- p %>% layout(
#   scene = list(
#     xaxis = list(title = "delta", showgrid = FALSE),
#     yaxis = list(title = "xi", showgrid = FALSE),
#     zaxis = list(title = "Density", showgrid = FALSE),
#     aspectratio = list(x = 1, y = 1, z = 0.7)
#   ),
#   title = "Joint Posterior Density"
# )
