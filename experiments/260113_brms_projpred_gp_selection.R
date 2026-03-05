library(brms)
library(projpred)

# dat must contain: data/simulated/simulated_data_P4_N50.R

# 1) Joint GP over (x1, x2, x3, x4): y = f(x1,x2,x3,x4) + error
fit_joint <- brm(
  y ~ 1 + gp(x1, x2, x3, x4),
  data = df,
  family = gaussian(),
  backend = "rstan",
  chains = 4,
  cores = 4,
  iter = 2000
)

# 2) Additive GP (sum of 1D GPs): y = f1(x1)+f2(x2)+f3(x3)+f4(x4) + error
fit_additive <- brm(
  y ~ 1 + gp(x1) + gp(x2) + gp(x3) + gp(x4),
  data = df,
  family = gaussian(),
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  iter = 2000
)

# Perform projpred variable selection on joint model
refm_joint <- get_refmodel(fit_joint)
