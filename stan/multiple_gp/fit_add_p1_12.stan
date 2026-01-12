data {
  int<lower=1> N;
  int<lower=1> P;                 // number of predictors
  matrix[N, P] X;                 // raw inputs (each column is a predictor)
  vector[N] y;                    // response (raw)
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
  // --- Per-predictor GP hyperparameters (delta on sigma-scale) ---
  vector<lower=0>[P] lambda;      // length-scales
  vector<lower=0>[P] delta;       // standardized GP amplitudes

  // Standardized linear coefficients: theta_j = beta1_j / sigma
  vector[P] theta;

  // noise
  real<lower=0> sigma;            // noise SD

  // intercept (centered scale)
  real mu;                        // mean of y_c
}

transformed parameters {
  vector<lower=0>[P] alpha = delta * sigma;

  // de-standardized slopes
  vector[P] beta1 = sigma * theta;

  // Prior buffer (for bridge sampling)
  real lprior = 0;

  // --- Priors on GP terms ---
  // delta_j ~ half-t(4, 2.7)
  for (j in 1:P)
    lprior += student_t_lpdf(delta[j] | 4, 0, 2.7) - student_t_lccdf(0 | 4, 0, 2.7);

  // lambda_j ~ IG(3, 1.5)
  for (j in 1:P)
    lprior += inv_gamma_lpdf(lambda[j] | 3, 1.5);

  // --- Shared priors ---
  // sigma ~ half-t(3, 2.5)
  lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);

  // mu ~ Student-t(3, 0, 10)
  lprior += student_t_lpdf(mu | 3, 0, 10);

  // theta_j ~ Cauchy(0, sqrt(N)/sqrt(x_j' x_j))   (JZS-style on standardized slopes)
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

  // ----- Accumulate orthogonalized GP covariances (NO linear marginalization) -----
  matrix[N, N] Sigma = rep_matrix(0, N, N);

  // Sum of Ktilde_j over predictors
  for (j in 1:P) {
    // Design for orthogonalization for predictor j: X_j = [1, x_cj]
    matrix[N, 2] Xj;
    Xj[, 1] = rep_vector(1.0, N);
    Xj[, 2] = X_c[, j];

    // Projector P_j onto span{1, x_cj}
    matrix[2, 2] XX = crossprod(Xj);
    matrix[2, N] S  = mdivide_left_spd(XX, Xj');
    matrix[N, N] Pj = Xj * S;

    matrix[N, N] I_minus_Pj = I - Pj;

    // Latent GP covariance for predictor j
    matrix[N, N] Kj = gp_exp_quad_cov(to_array_1d(X_c[, j]), alpha[j], lambda[j]);
    for (n in 1:N) Kj[n, n] += 1e-12;

    // Orthogonalized covariance: Ktilde_j = (I - Pj) Kj (I - Pj), symmetrized
    matrix[N, N] temp = I_minus_Pj * Kj * I_minus_Pj;
    matrix[N, N] Ktilde_j = 0.5 * (temp + temp');

    // Accumulate into Sigma
    Sigma += Ktilde_j;
  }

  // Add noise variance
  for (n in 1:N) Sigma[n, n] += square(sigma);

  // Mean on centered scale: intercept + linear part
  vector[N] mu_vec = rep_vector(mu, N) + X_c * beta1;

  // Likelihood
  {
    matrix[N, N] L = cholesky_decompose(Sigma);
    target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
  }
}

generated quantities {
  // Minimal artifacts for plotting partial effects on grids
  real           beta0;        // intercept on ORIGINAL scale (for grids)
  matrix[50, P]  x_grid;            // raw-scale grids per predictor
  matrix[50, P]  f_mean_grid;       // partial GP mean along x_j on the grid
  matrix[50, P]  linear_est_grid;   // beta0 + beta1_j * x_c_grid
  matrix[50, P]  full_est_grid;     // linear_est_grid + f_mean_grid

  {
    // ---------- Build Σ (GP only) and cache (I - Pj) ----------
    matrix[N, N] I = diag_matrix(rep_vector(1.0, N));
    matrix[N, N] Sigma = rep_matrix(0, N, N);
    array[P] matrix[N, N] IminPj_arr;

    for (j in 1:P) {
      // [1, x_cj] projector
      matrix[N, 2] Xj;
      Xj[, 1] = rep_vector(1.0, N);
      Xj[, 2] = X_c[, j];

      matrix[2, 2] XX = crossprod(Xj);
      matrix[2, N] S  = mdivide_left_spd(XX, Xj');
      matrix[N, N] Pj = Xj * S;
      matrix[N, N] I_minus_Pj = I - Pj;
      IminPj_arr[j] = I_minus_Pj;

      // GP covariance and orthogonalization
      array[N] real xj_arr = to_array_1d(X_c[, j]);
      matrix[N, N] Kj = gp_exp_quad_cov(xj_arr, delta[j] * sigma, lambda[j]);
      for (n in 1:N) Kj[n, n] += 1e-12;

      matrix[N, N] temp = I_minus_Pj * Kj * I_minus_Pj;
      matrix[N, N] Ktilde_j = 0.5 * (temp + temp');

      Sigma += Ktilde_j;
    }
    for (n in 1:N) Sigma[n, n] += square(sigma);

    // ---------- Precision action a = Σ^{-1} (y_c - mean) ----------
    matrix[N, N] L = cholesky_decompose(Sigma);
    vector[N] mean_vec = rep_vector(mu, N) + X_c * beta1;   // beta1 from transformed params
    vector[N] resid    = y_c - mean_vec;

    vector[N] v = mdivide_left_tri_low(L, resid);
    vector[N] a = mdivide_right_tri_low(v', L)';            // a = Σ^{-1} resid

    // ---------- Grids (partial effects only) ----------
    beta0 = mu + y_bar;

    for (j in 1:P) {
      // 50-point raw grid for x_j
      real xmin = min(col(X, j));
      real xmax = max(col(X, j));
      for (g in 1:50)
        x_grid[g, j] = xmin + (xmax - xmin) * (g - 1) / 49;

      // center the grid to match training centering
      vector[50] x_c_grid;
      for (g in 1:50) x_c_grid[g] = x_grid[g, j] - x_bar[j];

      // Grid design X_new = [1, x_c_grid]
      matrix[50, 2] X_new;
      X_new[, 1] = rep_vector(1.0, 50);
      X_new[, 2] = x_c_grid;

      // Grid projector and I - P for grid side
      matrix[50, 50] I_G = diag_matrix(rep_vector(1.0, 50));
      matrix[2, 2]   XX_new = crossprod(X_new);
      matrix[2, 2]   XX_inv = inverse_spd(XX_new);
      matrix[50, 50] P_new  = X_new * XX_inv * X_new';
      matrix[50, 50] IminP_n = I_G - P_new;

      // Training-side I - P for this predictor
      matrix[N, N] IminP_tr = IminPj_arr[j];

      // Cross-cov (grid × train) for component j, then orthogonalize both sides
      array[50] real x_c_grid_arr = to_array_1d(x_c_grid);
      array[N]  real xj_arr       = to_array_1d(X_c[, j]);
      matrix[50, N] K_cross = gp_exp_quad_cov(x_c_grid_arr, xj_arr,
                                              delta[j] * sigma, lambda[j]);
      matrix[50, N] Ktilde_cross = IminP_n * K_cross * IminP_tr;

      // Grid partial means and displays (ORIGINAL scale)
      f_mean_grid[, j]     = Ktilde_cross * a;
      linear_est_grid[, j] = rep_vector(beta0, 50) + beta1[j] * x_c_grid;
      full_est_grid[, j]   = linear_est_grid[, j] + f_mean_grid[, j];
    }
  }
}
