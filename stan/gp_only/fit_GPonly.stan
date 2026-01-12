functions {
  // 1D SE kernel for a single input vector (returns N x N)
  matrix gp_exp_quad_cov_vec(vector x, real alpha, real l) {
    int N = rows(x);
    array[N] real xa = to_array_1d(x);
    return gp_exp_quad_cov(xa, alpha, l);
  }
}

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

  // noise
  real<lower=0> sigma;            // noise SD

  // GP mean
  real mu;                        // mean of y_c
}

transformed parameters {
  vector<lower=0>[P] alpha = delta * sigma;

  // Prior buffer (for bridge sampling)
  real lprior = 0;

  // half-t priors via truncated Student-t
  // delta_j ~ half-t(4, 2.7)
  for (j in 1:P)
    lprior += student_t_lpdf(delta[j] | 4, 0, 2.7) - student_t_lccdf(0 | 4, 0, 2.7);

  // lambda_j ~ IG(3, 0.8)
  for (j in 1:P)
    lprior += inv_gamma_lpdf(lambda[j] | 3, 0.8);

  // sigma ~ half-t(3, 2.5)
  lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);

  // mu ~ Student-t(3, 0, 10)
  lprior += student_t_lpdf(mu | 3, 0, 10);
}

model {
  target += lprior;

  // ----- Sum of (non-orthogonal) GP covariances -----
  matrix[N, N] Sigma = rep_matrix(0, N, N);

  for (j in 1:P) {
    matrix[N, N] Kj = gp_exp_quad_cov_vec(X_c[, j], alpha[j], lambda[j]);
    for (n in 1:N) Kj[n, n] += 1e-12;   // jitter
    Sigma += Kj;
  }

  // Add noise variance
  for (n in 1:N) Sigma[n, n] += square(sigma);

  // Mean on centered scale: intercept only
  {
    matrix[N, N] L = cholesky_decompose(Sigma);
    vector[N] mu_vec = rep_vector(mu, N);
    target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
  }
}

generated quantities {
  // ----- Outputs on ORIGINAL scale -----
  matrix[N, P] f_mean_mat;      // training-point GP means per component
  vector[N]    f_mean_total;    // sum_j f_j at training points
  real         beta0_orig;      // intercept on ORIGINAL scale
  vector[N]    full_est;        // beta0_orig + sum_j f_j

  // Per-predictor 50-point grids (ORIGINAL scale)
  matrix[50, P] x_grid;         // raw-scale grids
  matrix[50, P] f_mean_grid;    // partial GP mean along x_j on the grid
  matrix[50, P] full_est_grid;  // beta0_orig + f_mean_grid

  {
    // ---------- Rebuild Σ and K_sum ----------
    matrix[N, N] Sigma = rep_matrix(0, N, N);
    matrix[N, N] K_sum = rep_matrix(0, N, N);

    for (j in 1:P) {
      matrix[N, N] Kj = gp_exp_quad_cov_vec(X_c[, j], alpha[j], lambda[j]);
      for (n in 1:N) Kj[n, n] += 1e-12;
      Sigma += Kj;
      K_sum += Kj;
    }
    for (n in 1:N) Sigma[n, n] += square(sigma);

    // a = Σ^{-1} (y_c - mu)
    matrix[N, N] L = cholesky_decompose(Sigma);
    vector[N] resid = y_c - rep_vector(mu, N);
    vector[N] v = mdivide_left_tri_low(L, resid);
    vector[N] a = mdivide_right_tri_low(v', L)';   // a = Σ^{-1} resid

    // ---------- Training-point component means ----------
    for (j in 1:P) {
      matrix[N, N] Kj = gp_exp_quad_cov_vec(X_c[, j], alpha[j], lambda[j]);
      for (n in 1:N) Kj[n, n] += 1e-12;
      f_mean_mat[, j] = Kj * a;                 // Cov(f_j, y) Σ^{-1} resid
    }

    // totals and fitted values (ORIGINAL scale)
    f_mean_total = rep_vector(0.0, N);
    for (j in 1:P) f_mean_total += f_mean_mat[, j];

    beta0_orig = mu + y_bar;
    full_est   = beta0_orig + f_mean_total;

    // ---------- Grids (partial effects only) ----------
    for (j in 1:P) {
      // 50-point raw grid for x_j
      real xmin = min(col(X, j));
      real xmax = max(col(X, j));
      for (g in 1:50)
        x_grid[g, j] = xmin + (xmax - xmin) * (g - 1) / 49;

      // center grid to match training centering
      vector[50] x_c_grid;
      for (g in 1:50) x_c_grid[g] = x_grid[g, j] - x_bar[j];

      // Cross-cov (grid × train), NO orthogonalization
      array[50] real x_c_grid_arr = to_array_1d(x_c_grid);
      array[N]  real xj_arr       = to_array_1d(X_c[, j]);
      matrix[50, N] K_cross = gp_exp_quad_cov(x_c_grid_arr, xj_arr,
                                              delta[j] * sigma, lambda[j]);

      // Partial GP mean and full grid curve (ORIGINAL scale)
      f_mean_grid[, j]   = K_cross * a;
      full_est_grid[, j] = rep_vector(beta0_orig, 50) + f_mean_grid[, j];
    }
  }
}
