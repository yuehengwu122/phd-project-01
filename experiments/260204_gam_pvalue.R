# Load required libraries
library(gam) # For Generalized Additive Models (GAM)
library(mgcv) # For modern GAM fitting with penalized splines

# -----------------------------------------------------------------------------
# Simulation Study for Type I Error Rate
# -----------------------------------------------------------------------------

library(mgcv)
set.seed(123)
n <- 50
beta <- 0.5

# Generate p-values from GAM fits
pvals <- unlist(lapply(1:1e3, function(s) {
  x <- seq(-2, 2, length.out = n)
  y <- x * beta + rnorm(n)
  gam1 <- mgcv::gam(y ~ x + s(x))
  sum1 <- summary(gam1)
  # Get p-value of nonlinear effect
  sum1$s.pv
}))

# Filter out NA and non-positive p-values
pvals <- pvals[!is.na(pvals) & pvals > 0]

# Plot histogram of p-values
hist(pvals, freq = FALSE, main = "Histogram of p-values", xlab = "p-value")

# Calculate inflated type I error rate at significance level 0.05
mean(pvals < 0.05)


pvals_penalty <- unlist(lapply(1:1e3, function(s) {
  x <- seq(-2, 2, length.out = n)
  y <- x * beta + rnorm(n)
  gam1 <- mgcv::gam(y ~ x + s(x, bs = 'tp', m = 1, k = 10))
  sum1 <- summary(gam1)
  # Get p-value of nonlinear effect
  sum1$s.pv
}))

pvals_penalty <- pvals_penalty[!is.na(pvals_penalty) & pvals_penalty > 0]

hist(
  pvals_penalty,
  freq = FALSE,
  main = "Histogram of p-values with penalty",
  xlab = "p-value"
)

mean(pvals_penalty < 0.05)


pvals_gam <- unlist(lapply(1:1e3, function(s) {
  x <- seq(-2, 2, length.out = n)
  y <- x * beta + rnorm(n)
  gam1 <- gam::gam(y ~ x + s(x))
  sum1 <- summary(gam1)
  # Get p-value of nonlinear effect
  sum1$anova$`Pr(F)`[3]
}))

pvals_gam <- pvals_gam[!is.na(pvals_gam) & pvals_gam > 0]
hist(
  pvals_gam,
  freq = FALSE,
  main = "Histogram of p-values (gam package)",
  xlab = "p-value"
)
mean(pvals_gam < 0.05)


# -----------------------------------------------------------------------------
# Compare p-values for linear and nonlinear effects between `gam` and `mgcv` packages
# -----------------------------------------------------------------------------

# Generate p-values for linear and nonlinear effects using `mgcv`
pvals_mgcv <- lapply(1:1e3, function(s) {
  x <- seq(-2, 2, length.out = n)
  y <- sin(2 * x) + rnorm(n) # Strong nonlinear effect of x on y, no linear effect
  gam1 <- mgcv::gam(y ~ x + s(x, bs = 'tp', m = 1, k = 10))
  sum1 <- summary(gam1)
  list(
    linear_p = sum1$p.table["x", "Pr(>|t|)"],
    nonlinear_p = sum1$s.table["s(x)", "p-value"]
  )
})

# Extract and filter p-values for linear and nonlinear effects
pvals_linear_mgcv <- unlist(lapply(pvals_mgcv, `[[`, "linear_p"))
pvals_linear_mgcv <- pvals_linear_mgcv[
  !is.na(pvals_linear_mgcv) & pvals_linear_mgcv > 0
]

pvals_nonlinear_mgcv <- unlist(lapply(pvals_mgcv, `[[`, "nonlinear_p"))
pvals_nonlinear_mgcv <- pvals_nonlinear_mgcv[
  !is.na(pvals_nonlinear_mgcv) & pvals_nonlinear_mgcv > 0
]

# Generate p-values for linear and nonlinear effects using `gam`
pvals_gam <- lapply(1:1e3, function(s) {
  x <- seq(-2, 2, length.out = n)
  y <- sin(2 * x) + rnorm(n) # Strong nonlinear effect of x on y, no linear effect
  gam1 <- gam::gam(y ~ x + s(x))
  sum1 <- summary(gam1)
  list(
    linear_p = sum1$parametric.anova$`Pr(>F)`[1],
    nonlinear_p = sum1$anova$`Pr(F)`[3]
  )
})

# Extract and filter p-values for linear and nonlinear effects
pvals_linear_gam <- unlist(lapply(pvals_gam, `[[`, "linear_p"))
pvals_linear_gam <- pvals_linear_gam[
  !is.na(pvals_linear_gam) & pvals_linear_gam > 0
]

pvals_nonlinear_gam <- unlist(lapply(pvals_gam, `[[`, "nonlinear_p"))
pvals_nonlinear_gam <- pvals_nonlinear_gam[
  !is.na(pvals_nonlinear_gam) & pvals_nonlinear_gam > 0
]

# Compare histograms of p-values
par(mfrow = c(1, 2))
hist(
  pvals_linear_mgcv,
  freq = FALSE,
  main = "mgcv: Linear Effect (Null)",
  xlab = "p-value"
)
hist(
  pvals_linear_gam,
  freq = FALSE,
  main = "gam: Linear Effect (Null)",
  xlab = "p-value"
)

# Compare Type I error rates at significance level 0.05
cat("Type I error rate (mgcv):", mean(pvals_linear_mgcv < 0.05), "\n")
cat("Type I error rate (gam):", mean(pvals_linear_gam < 0.05), "\n")

# Compare power to detect nonlinear effect at significance level 0.05
cat(
  "Power to detect nonlinear effect (mgcv):",
  mean(pvals_nonlinear_mgcv < 0.05),
  "\n"
)
cat(
  "Power to detect nonlinear effect (gam):",
  mean(pvals_nonlinear_gam < 0.05),
  "\n"
)


# -----------------------------------------------------------------------------
# General Function for Bootstrap Analysis of p-values
# -----------------------------------------------------------------------------

bootstrap_gam_analysis <- function(gam_data, dataset_name, n_bootstrap = 100) {
  set.seed(123)
  predictors <- colnames(gam_data)[-1] # Get predictor names

  # Initialize results storage as matrices with zeros
  results_mgcv <- matrix(0, nrow = length(predictors), ncol = 2)
  results_gam <- matrix(0, nrow = length(predictors), ncol = 2)
  colnames(results_mgcv) <- c("linear_sig_rate", "nonlinear_sig_rate")
  colnames(results_gam) <- c("linear_sig_rate", "nonlinear_sig_rate")
  rownames(results_mgcv) <- predictors
  rownames(results_gam) <- predictors

  cat(
    "Starting bootstrap analysis for",
    dataset_name,
    "with",
    n_bootstrap,
    "replicates...\n"
  )
  cat("This may take a few minutes...\n\n")

  # Create full model formulas
  formula_mgcv <- as.formula(paste0(
    "y ~ ",
    paste(predictors, collapse = " + "),
    " + ",
    paste0("s(", predictors, ", bs='tp', m=1, k=10)", collapse = " + ")
  ))

  formula_gam <- as.formula(paste0(
    "y ~ ",
    paste(predictors, collapse = " + "),
    " + ",
    paste0("s(", predictors, ")", collapse = " + ")
  ))

  # Bootstrap analysis with full models
  for (b in 1:n_bootstrap) {
    cat("Bootstrap replicate:", b, "\n")

    # Bootstrap sampling
    boot_indices <- sample(nrow(gam_data), replace = TRUE)
    boot_data <- gam_data[boot_indices, ]

    # mgcv fit
    tryCatch(
      {
        gam_mgcv <- mgcv::gam(formula_mgcv, data = boot_data)
        sum_mgcv <- summary(gam_mgcv)

        # Extract p-values for each predictor
        for (pred in predictors) {
          linear_p <- sum_mgcv$p.table[pred, "Pr(>|t|)"]
          nonlinear_p <- sum_mgcv$s.table[paste0("s(", pred, ")"), "p-value"]

          if (!is.na(linear_p) && linear_p > 0) {
            results_mgcv[pred, "linear_sig_rate"] <- results_mgcv[
              pred,
              "linear_sig_rate"
            ] +
              (linear_p < 0.05)
          }
          if (!is.na(nonlinear_p) && nonlinear_p > 0) {
            results_mgcv[pred, "nonlinear_sig_rate"] <- results_mgcv[
              pred,
              "nonlinear_sig_rate"
            ] +
              (nonlinear_p < 0.05)
          }
        }
      },
      error = function(e) {
        cat("mgcv fit failed for replicate", b, "\n")
      }
    )

    # gam fit
    tryCatch(
      {
        gam_fit <- gam::gam(formula_gam, data = boot_data)
        sum_gam <- summary(gam_fit)

        # Extract p-values for each predictor
        for (i in seq_along(predictors)) {
          pred <- predictors[i]
          linear_p <- sum_gam$parametric.anova$`Pr(>F)`[i]
          nonlinear_p <- sum_gam$anova$`Pr(F)`[length(predictors) + i + 1] # Adjust index for smooth terms

          if (!is.na(linear_p) && linear_p > 0) {
            results_gam[pred, "linear_sig_rate"] <- results_gam[
              pred,
              "linear_sig_rate"
            ] +
              (linear_p < 0.05)
          }
          if (!is.na(nonlinear_p) && nonlinear_p > 0) {
            results_gam[pred, "nonlinear_sig_rate"] <- results_gam[
              pred,
              "nonlinear_sig_rate"
            ] +
              (nonlinear_p < 0.05)
          }
        }
      },
      error = function(e) {
        cat("gam fit failed for replicate", b, "\n")
      }
    )
  }

  # Convert counts to proportions
  results_mgcv <- results_mgcv / n_bootstrap
  results_gam <- results_gam / n_bootstrap

  # Print results
  cat("\n=== ", dataset_name, " Results Summary ===\n\n")
  cat("mgcv package results:\n")
  print(round(results_mgcv, 3))
  cat("\ngam package results:\n")
  print(round(results_gam, 3))

  return(list(mgcv = results_mgcv, gam = results_gam))
}

# -----------------------------------------------------------------------------
# Apply Analysis to Different Datasets
# -----------------------------------------------------------------------------

# Diabetes Dataset
data(efron2004, package = "care")
x_mat <- as.matrix(efron2004$x)
class(x_mat) <- "matrix"
X <- x_mat[, -2] # Remove "sex" (a binary variable)

diabetes_data <- data.frame(
  y = efron2004$y[, 1],
  x_mat[, -2]
)
colnames(diabetes_data)[-1] <- colnames(efron2004$x)[-2]

# Run analysis for diabetes dataset
diabetes_results <- bootstrap_gam_analysis(diabetes_data, "Diabetes Dataset")

# Boston Housing Dataset
data("Boston", package = "MASS")
boston_data <- data.frame(
  y = Boston$medv, # Median value of owner-occupied homes
  as.matrix(Boston[, -c(4, 9, 14)]) # Exclude categorical predictors
)

# Run analysis for Boston Housing dataset
boston_results <- bootstrap_gam_analysis(boston_data, "Boston Housing Dataset")

lm(y ~ age + bmi + bp + s1 + s2 + s3 + s4 + s5 + s6, data = gam_data)
