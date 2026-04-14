data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;
  vector[N] y;
}

transformed data {
  // ----- Center y and each column of X -----
  real y_bar = mean(y);
  vector[N] y_c = y - y_bar;

  row_vector[P] x_bar;
  matrix[N, P] X_c;
  for (j in 1:P) {
    x_bar[j] = mean(col(X, j));
    X_c[, j] = col(X, j) - x_bar[j];
  }

  // ---- objects for mixture-of-g prior on theta ----
  matrix[P, P] XtX = crossprod(X_c);
  matrix[P, P] XtX_ridge = XtX + 1e-10 * diag_matrix(rep_vector(1.0, P));
  matrix[P, P] V_theta_base = inverse_spd(XtX_ridge);
  matrix[P, P] L_V_theta_base = cholesky_decompose(V_theta_base);
}

parameters {
  real log_g_theta;
  vector[P] z_theta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[P] theta = exp(0.5 * log_g_theta) * (L_V_theta_base * z_theta);
  vector[P] beta1 = sigma * theta;

  real lprior = 0;

  // sigma prior
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
            - student_t_lccdf(0 | 3, 0, 2.5);

  // log_g_theta prior
  lprior += normal_lpdf(log_g_theta | log(N), 0.5);
}

model {
  target += lprior;
  z_theta ~ normal(0, 1);

  // Sigma = sigma^2 * I
  real s2 = square(sigma);
  matrix[N, N] Sigma = diag_matrix(rep_vector(s2, N));

  // Likelihood
  vector[N] mu_vec = X_c * beta1;
  {
    matrix[N, N] L = cholesky_decompose(Sigma);
    target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
  }
}
