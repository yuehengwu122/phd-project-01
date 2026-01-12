# ============================================================
# Compile Stan Models
# ============================================================

# Main model

gp_add_p1_11 <- rstan::stan_model(file = "stan/multiple_gp/fit_add_p1_11.stan")

gp_add_p1_12 <- rstan::stan_model(file = "stan/multiple_gp/fit_add_p1_12.stan")

gp_add_p1_13 <- rstan::stan_model(file = "stan/multiple_gp/fit_add_p1_13.stan")

gp_add_p1_14 <- rstan::stan_model(file = "stan/multiple_gp/fit_add_p1_14.stan")

gp_add_p1_15 <- rstan::stan_model(file = "stan/multiple_gp/fit_add_p1_15.stan")


# Null model for Bayes Factor
gp_null <- rstan::stan_model(file = "stan/multiple_gp/fit_null.stan")

# single GP model in simulation study
gp_single_p1 <- rstan::stan_model(file = "stan/single_gp/fit_single_07.stan")
gp_single_p2 <- rstan::stan_model(file = "stan/single_gp/fit_single_07_p2.stan")
gp_single_null <- rstan::stan_model(file = "stan/single_gp/null_4.stan")