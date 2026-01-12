data {
  int<lower=1> N;
  int<lower=1> P; // number of predictors
  matrix[N, P] X; // raw inputs (each column is a predictor)
  vector[N] y; // response (raw)
}

transformed data {
  // ----- Center y and each column of X -----
  real y_bar = mean(y);
  vector[N] y_c = y - y_bar;

  // column means of X and centered matrix
  row_vector[P] x_bar;
  matrix[N, P] X_c;
  for (j in 1:P) {
    x_bar[j] = mean(col(X, j));
    X_c[, j] = col(X, j) - x_bar[j];
  }
}

parameters {
  // Standardized linear coefficients: theta_j = beta_j / sigma
  vector[P] theta;

  // noise
  real<lower=0> sigma; // noise SD

  // intercept (centered scale)
  real mu; // mean of y_c
}

transformed parameters {
  // de-standardized slopes
  vector[P] beta = sigma * theta;

  // Prior buffer (for bridge sampling)
  real lprior = 0;

  // --- Shared priors ---
  // sigma ~ half-t(3, 2.5)
  lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);

  // mu ~ Student-t(3, 0, 10)
  lprior += student_t_lpdf(mu | 3, 0, 10);

  // theta_j ~ Cauchy(0, sqrt(N)/sqrt(x_j' x_j))   (JZS-style on standardized slopes)
  for (j in 1:P) {
    real xjxj = dot_self(X_c[,j]); // x_j'x_j
    real scale_theta = sqrt(N) / sqrt(xjxj);
    lprior += cauchy_lpdf(theta[j] | 0.0, scale_theta);
  }
}

model {
  target += lprior;

  // ----- Common identity -----
  matrix[N, N] I = diag_matrix(rep_vector(1.0, N));

  matrix[N, N] Sigma = rep_matrix(0, N, N);

  // Add noise variance
  for (n in 1:N) Sigma[n, n] += square(sigma);

  // Mean on centered scale: intercept + linear part
  vector[N] mu_vec = rep_vector(mu, N) + X_c * beta;

  // Likelihood
  {
    matrix[N, N] L = cholesky_decompose(Sigma);
    target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
  }
}
