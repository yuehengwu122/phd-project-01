source("_common.R")

# ============================================================
# Part 1: Simulate Data with Known GP Effect
# ============================================================
set.seed(123)
N <- 50
X <- seq(-0.5, 0.5, length = N) |> matrix(ncol = 1)

# Generate data with
# (1) true delta = 3 (nonlinear effect exists)
# (2) true delta = 0.1 (no nonlinear effect)
y <- generate_gp_data(
  delta = 3,
  lambda = 0.3,
  sigma = 0.1,
  N = 50,
  x = X,
  n_samples = 1
) |>
  as.vector()

# ============================================================
# Part 2: Fit GP Models with Two Different Priors
# ============================================================

# Fit with p1 prior (delta ~ half-t(4, 2.7))
fit_gp_p1 <- fit_model(
  gp_add_p1_15,
  X = X,
  y = y,
  chains = 2,
  iter = 5e3,
  warmup = 500
)

fit_gp_p2 <- fit_model(
  gp_add_p2_15,
  X = X,
  y = y,
  chains = 2,
  iter = 5e3,
  warmup = 500
)

fit_gp_p3 <- fit_model(
  gp_add_p3_15,
  X = X,
  y = y,
  chains = 2,
  iter = 5e3,
  warmup = 500
)

fit_gp_p4 <- fit_model(
  gp_add_p4_15,
  X = X,
  y = y,
  chains = 2,
  iter = 5e3,
  warmup = 500
)

fit_null <- fit_model(
  gp_null,
  X = X,
  y = y,
  chains = 2,
  iter = 5e3,
  warmup = 500
)

bf_results <- tibble::tibble(
  prior = c("p1", "p2", "p3", "p4"),
  bayes_factor = c(
    compute_bayes_factor_bridge(fit_gp_p1, fit_null, log = TRUE),
    compute_bayes_factor_bridge(fit_gp_p2, fit_null, log = TRUE),
    compute_bayes_factor_bridge(fit_gp_p3, fit_null, log = TRUE),
    compute_bayes_factor_bridge(fit_gp_p4, fit_null, log = TRUE)
  )
)

# Print results as a table with 4 digits
bf_results |>
  dplyr::mutate(bayes_factor = round(bayes_factor, 4)) |>
  print()

post_delta_k <- data.frame(
  p1 = as.vector(rstan::extract(fit_gp_p1$fit)$delta),
  p2 = as.vector(rstan::extract(fit_gp_p2$fit)$delta2),
  p3 = as.vector(rstan::extract(fit_gp_p3$fit)$delta3),
  p4 = as.vector(rstan::extract(fit_gp_p4$fit)$delta4)
)

# Compute BF10s for various eps and d values
lapply(1:4, function(i) {
  bf10_interval_null_delta_1d(
    post_delta_k[[i]],
    prior_id = paste0("p", i),
    eps = c(0, 0.1, 0.3, 0.5, 0.7, 1)^i,
    d = c(2, 4, 6)
  )
})


boxplot(
  post_delta_k,
  names = c("p1", "p2", "p3", "p4"),
  main = "Boxplot of Delta^k Values for Different Priors",
  ylab = "Delta",
  ylim = c(0, 100)
)

boxplot(
  post_delta,
  names = c("p1", "p2", "p3", "p4"),
  main = "Boxplot of Delta Values for Different Priors",
  ylab = "Delta"
)


compute_sdr_from_fit(fit_a3_dep, prior_id = "p1")
compute_sdr_from_fit(fit_a3_dep, prior_id = "p1", d = 4)

compute_sdr_from_fit(fit_a4_dep, prior_id = "p1")
compute_sdr_from_fit(fit_a4_dep, prior_id = "p1", d = 4)







library(dplyr)

# delta-scale grid
eps_grid <- seq(0, 1, length.out = 200)

# compute across k=1..4
res_grid <- bind_rows(lapply(1:4, function(k) {

  out <- bf10_interval_null_delta_1d(
    post_delta_k[[k]],
    prior_id = paste0("p", k),
    eps = eps_grid^k
  )

  out %>%
    mutate(
      k = k,
      eps_delta = epsilon^(1 / k), # back to delta-scale
      log_bf10 = log_bf10 # bf10 already on log scale
    )
}))

# base R plot
yl <- range(res_grid$log_bf10, finite = TRUE)

plot(NA,
     xlim = c(0, 1),
     ylim = yl,
     xlab = expression(epsilon~"(delta-scale)"),
     ylab = "log BF10",
     main = "Interval-null log BF10 vs epsilon")

for (k in 1:4) {
  tmp <- res_grid %>% filter(k == !!k)
  lines(tmp$eps_delta, tmp$log_bf10, lwd = 2, col = k)
}

legend("topright",
       legend = paste0("k = ", 1:4),
       col = 1:4,
       lwd = 2,
       bty = "n")
