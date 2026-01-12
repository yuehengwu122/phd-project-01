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
  
  matrix[N,N] I  = diag_matrix(rep_vector(1.0,N));
  vector[N]  one = rep_vector(1.0, N);
  matrix[N,N] P0 = (one * one') / N;       // projector onto constants
  matrix[N,N] I0 = I - P0;

  for (j in 1:P) {
    matrix[N,N] Kj = gp_exp_quad_cov_vec(X_c[, j], alpha[j], lambda[j]);
    for (n in 1:N) Kj[n,n] += 1e-12;
    matrix[N,N] Ktilde_j = 0.5 * (I0 * Kj * I0 + (I0 * Kj * I0)');  // symmetrize
    Sigma += Ktilde_j;
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
  real           beta0_orig;        // intercept on ORIGINAL scale
  matrix[50, P]  x_grid;            // raw-scale grids per predictor
  matrix[50, P]  f_mean_grid;       // partial GP mean along x_j (constant-orthogonalized)
  matrix[50, P]  full_est_grid;     // beta0_orig + f_mean_grid

  {
    // ----- Constant projector on training set -----
    matrix[N, N] I   = diag_matrix(rep_vector(1.0, N));
    vector[N]    one = rep_vector(1.0, N);
    matrix[N, N] P0  = (one * one') / N;   // projector onto constants
    matrix[N, N] I0  = I - P0;             // remove constant direction

    // ----- Build Σ = sum_j Ktilde_j + σ^2 I -----
    matrix[N, N] Sigma = rep_matrix(0, N, N);
    for (j in 1:P) {
      array[N]  real xj_arr = to_array_1d(X_c[, j]);
      matrix[N, N] Kj = gp_exp_quad_cov(xj_arr, delta[j] * sigma, lambda[j]);
      for (n in 1:N) Kj[n, n] += 1e-12;

      matrix[N, N] temp = I0 * Kj * I0;
      matrix[N, N] Ktilde_j = 0.5 * (temp + temp');  // symmetrize
      Sigma += Ktilde_j;
    }
    for (n in 1:N) Sigma[n, n] += square(sigma);

    // ----- a = Σ^{-1} (y_c - mu) -----
    matrix[N, N] L = cholesky_decompose(Sigma);
    vector[N] resid = y_c - rep_vector(mu, N);
    vector[N] v = mdivide_left_tri_low(L, resid);
    vector[N] a = mdivide_right_tri_low(v', L)';   // a = Σ^{-1} resid

    // ----- Intercept on ORIGINAL scale -----
    beta0_orig = mu + y_bar;

    // ----- Per-predictor grids (with constant removal on both sides) -----
    for (j in 1:P) {
      // 50-point raw grid for x_j
      real xmin = min(col(X, j));
      real xmax = max(col(X, j));
      for (g in 1:50)
        x_grid[g, j] = xmin + (xmax - xmin) * (g - 1) / 49;

      // center to match training centering
      vector[50] x_c_grid;
      for (g in 1:50) x_c_grid[g] = x_grid[g, j] - x_bar[j];

      // Grid-side constant projector
      matrix[50, 50] I_G  = diag_matrix(rep_vector(1.0, 50));
      vector[50]     oneG = rep_vector(1.0, 50);
      matrix[50, 50] P0G  = (oneG * oneG') / 50;
      matrix[50, 50] I0G  = I_G - P0G;

      // Cross-covariance K(grid, train), then remove constants on both sides
      array[50] real x_c_grid_arr = to_array_1d(x_c_grid);
      array[N]  real xj_arr       = to_array_1d(X_c[, j]);
      matrix[50, N] K_cross = gp_exp_quad_cov(x_c_grid_arr, xj_arr,
                                              delta[j] * sigma, lambda[j]);
      matrix[50, N] Ktilde_cross = I0G * K_cross * I0;

      // Partial GP mean (constant-orthogonalized) and full grid curve
      f_mean_grid[, j]   = Ktilde_cross * a;
      full_est_grid[, j] = rep_vector(beta0_orig, 50) + f_mean_grid[, j];
    }
  }
}
