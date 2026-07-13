# ============================================================
# Application 4: Boston Housing
#
# ARD baseline via gpboost
# Dataset: MASS::Boston
# - Predictors: 11 continuous features (excluding chas, rad)
# - Response:   medv (median home value)
# ============================================================

if (!requireNamespace("gpboost", quietly = TRUE)) {
  stop("Package 'gpboost' is required. Install it with install.packages('gpboost').")
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.")
}

data(Boston, package = "MASS")

X_raw <- as.matrix(Boston[, -c(4, 9, 14), drop = FALSE])  # drop chas, rad, medv
y <- Boston$medv

cat("N =", nrow(X_raw), ", P =", ncol(X_raw), "\n")
cat("Predictors:", paste(colnames(X_raw), collapse = ", "), "\n")

X_scaled <- scale(X_raw)
X_scaled <- as.matrix(X_scaled)

fit_gpboost <- gpboost::fitGPModel(
  y = y,
  X = matrix(1, nrow = nrow(X_scaled), ncol = 1),
  gp_coords = X_scaled,
  cov_function = "gaussian_ard",
  likelihood = "gaussian"
)

cov_pars <- gpboost::get_cov_pars(fit_gpboost)
range_idx <- grep("^GP_range_", names(cov_pars))
range_ard <- unname(cov_pars[range_idx])

ard_table <- data.frame(
  predictor = colnames(X_raw),
  range_ard = range_ard,
  inv_range_ard = 1 / range_ard,
  stringsAsFactors = FALSE
)

ard_by_range <- ard_table[order(ard_table$range_ard), ]
ard_by_relevance <- ard_table[order(-ard_table$inv_range_ard), ]

cat("\nARD ranges (ascending range)\n")
print(ard_by_range, row.names = FALSE)

cat("\nARD relevance ranking (descending 1 / range)\n")
print(ard_by_relevance, row.names = FALSE)

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

saveRDS(fit_gpboost, "results/fit_04_boston_gpboost_ard.rds")
write.csv(ard_table, "results/tab_04_boston_gpboost_ard.csv", row.names = FALSE)
