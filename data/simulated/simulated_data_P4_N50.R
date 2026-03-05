# Simulated Data Generation

source("R/data_prep.R")

# Set seed for reproducibility
set.seed(123)
n <- 50 # sample size

# 1. Simulate four predictors
x1 <- runif(n, -3, 3) # predictor with strong nonlinear effect
x2 <- runif(n, -3, 3) # predictor with moderate nonlinear effect
x3 <- runif(n, -3, 3) # predictor with linear effect
x4 <- runif(n, -3, 3) # predictor with no effect

# 2. Define effects on y
# Strong nonlinear: e.g. y += sin(π * x1^3)
# Moderate nonlinear: y += 0.5 * x2^2
# Linear effect: coef = 1.5
# No effect: x4 not used
y <- sin(2 * x1) *
  5 +
  -2 * x2 +
  0.8 * x2^2 +
  1.5 * x3 +
  rnorm(n, mean = 0, sd = 0.1) # noise

# Multivariate model
X <- as.matrix(cbind(
  rescale_to_range1(x1),
  rescale_to_range1(x2),
  rescale_to_range1(x3),
  rescale_to_range1(x4)
))

fit_data <- list(N = 50, y = y, X = X, P = 4)


# -----------------------#
# ---- Visualization ----#
# -----------------------#
library(ggplot2)
library(gridExtra)

df <- data.frame(x1, x2, x3, x4, y)

p1 <- ggplot(df, aes(x = x1, y = y)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "gam") +
  labs(title = "x1: strong nonlinear")

p2 <- ggplot(df, aes(x = x2, y = y)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "gam") +
  labs(title = "x2: linear + moderate nonlinear")

p3 <- ggplot(df, aes(x = x3, y = y)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "gam") +
  labs(title = "x3: linear")

p4 <- ggplot(df, aes(x = x4, y = y)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "gam") +
  labs(title = "x4: no effect")

grid.arrange(p1, p2, p3, p4, ncol = 2)
