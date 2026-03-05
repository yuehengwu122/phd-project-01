library(brms)

data("Boston", package = "MASS")
df <- Boston[, c('rm', 'medv')]
plot(df$rm, df$medv, pch = 16)


fit_gp_1 <- brm(
  medv ~ 1 + gp(rm),
  data = df,
  family = gaussian(),
  chains = 2,
  cores = 4
)

sc <- brms::make_stancode(
  medv ~ 1 + rm + gp(rm),
  data = df,
  family = gaussian()
)
cat(sc)
