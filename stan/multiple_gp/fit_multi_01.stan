functions {
  // 1D SE kernel for a single input vector (kept for parity; not used in ARD build)
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

  // ----- Global projector onto span{1, X_c} -----
  matrix[N, P + 1] Z;
  Z[, 1] = rep_vector(1.0, N);
  Z[, 2:(P + 1)] = X_c;

  matrix[P + 1, P + 1] ZZ = crossprod(Z);
  matrix[P + 1, N]     S  = mdivide_left_spd(ZZ, Z');
  matrix[N, N]         P_lin = Z * S;

  matrix[N, N] I = diag_matrix(rep_vector(1.0, N));
  matrix[N, N] IminP_lin = I - P_lin;
}

parameters {
  // --- ARD length-scales per predictor ---
  vector<lower=0>[P] lambda;      // length-scales (ARD)

  // --- Single GP amplitude on sigma-scale ---
  real<lower=0> delta;            // shared standardized GP amplitude

  // Standardized linear coefficients: theta_j = beta1_j / sigma
  vector[P] theta;

  // noise
  real<lower=0> sigma;            // noise SD

  // intercept (centered scale)
  real mu;                        // mean of y_c
}

transformed parameters {
  real<lower=0> alpha = delta * sigma;  // shared amplitude

  // de-standardized slopes
  vector[P] beta1 = sigma * theta;

  // Prior buffer (for bridge sampling)
  real lprior = 0;

  // --- Priors on GP terms ---
  // delta ~ half-t(4, 2.7)
  lprior += student_t_lpdf(delta | 4, 0, 2.7) - student_t_lccdf(0 | 4, 0, 2.7);

  // lambda_j ~ IG(3, 0.8)
  for (j in 1:P)
    lprior += inv_gamma_lpdf(lambda[j] | 3, 0.8);

  // --- Shared priors ---
  // sigma ~ half-t(3, 2.5)
  lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);

  // mu ~ Student-t(3, 0, 10)
  lprior += student_t_lpdf(mu | 3, 0, 10);

  // theta_j ~ Cauchy(0, 1/sqrt(x_j' x_j))   (JZS-style on standardized slopes)
  for (j in 1:P) {
    real scale_theta = 1.0 / sqrt(dot_self(X_c[, j]));
    lprior += cauchy_lpdf(theta[j] | 0.0, scale_theta);
  }
}

model {
  target += lprior;

  // ----- Build single ARD GP covariance K (N x N) -----
  matrix[N, N] K;
  {
    real a2 = square(alpha);
    for (i in 1:N) {
      // diagonal
      K[i, i] = a2 + 1e-12;
      for (n in (i + 1):N) {
        real r2 = 0;
        for (j in 1:P) {
          real diff = X_c[i, j] - X_c[n, j];
          r2 += (diff * diff) / square(lambda[j]);
        }
        real k = a2 * exp(-0.5 * r2);
        K[i, n] = k;
        K[n, i] = k;
      }
    }
  }

  // Orthogonalize to span{1, X_c}: Ktilde = (I - P_lin) K (I - P_lin), symmetrized
  matrix[N, N] temp = IminP_lin * K * IminP_lin;
  matrix[N, N] Ktilde = 0.5 * (temp + temp');

  // ----- Observation covariance Σ = Ktilde + σ^2 I -----
  matrix[N, N] Sigma = Ktilde;
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
  // Minimal artifacts for plotting partial effects on grids (partial dependence)
  real           beta0_orig;        // intercept on ORIGINAL scale (for grids)
  matrix[50, P]  x_grid;            // raw-scale grids per predictor
  matrix[50, P]  f_mean_grid;       // PD GP mean along x_j (others fixed at 0)
  matrix[50, P]  linear_est_grid;   // beta0_orig + beta1_j * x_c_grid
  matrix[50, P]  full_est_grid;     // linear_est_grid + f_mean_grid

  {
    // ---------- Rebuild ARD K and Σ (same as in model) ----------
    matrix[N, N] K;
    {
      real a2 = square(alpha);
      for (i in 1:N) {
        K[i, i] = a2 + 1e-12;
        for (n in (i + 1):N) {
          real r2 = 0;
          for (j in 1:P) {
            real diff = X_c[i, j] - X_c[n, j];
            r2 += (diff * diff) / square(lambda[j]);
          }
          real k = a2 * exp(-0.5 * r2);
          K[i, n] = k;
          K[n, i] = k;
        }
      }
    }
    matrix[N, N] temp = IminP_lin * K * IminP_lin;
    matrix[N, N] Ktilde = 0.5 * (temp + temp');

    matrix[N, N] Sigma = Ktilde;
    for (n in 1:N) Sigma[n, n] += square(sigma);

    // ---------- Precision action a = Σ^{-1} (y_c - mean) ----------
    matrix[N, N] L = cholesky_decompose(Sigma);
    vector[N] mean_vec = rep_vector(mu, N) + X_c * (sigma * theta);  // beta1
    vector[N] resid    = y_c - mean_vec;

    vector[N] v = mdivide_left_tri_low(L, resid);
    vector[N] a = mdivide_right_tri_low(v', L)';                      // a = Σ^{-1} resid

    // ---------- Grids: partial dependence by varying x_j, others fixed at 0 ----------
    beta0_orig = mu + y_bar;

    for (j in 1:P) {
      // 50-point raw grid for x_j
      real xmin = min(col(X, j));
      real xmax = max(col(X, j));
      for (g in 1:50)
        x_grid[g, j] = xmin + (xmax - xmin) * (g - 1) / 49;

      // center the grid to match training centering
      vector[50] x_c_grid;
      for (g in 1:50) x_c_grid[g] = x_grid[g, j] - x_bar[j];

      // Grid-side projector removing [1, x_c_grid] so display is pure deviation
      matrix[50, 2] X_new;
      X_new[, 1] = rep_vector(1.0, 50);
      X_new[, 2] = x_c_grid;

      matrix[2, 2]   XX_new  = crossprod(X_new);
      matrix[2, 2]   XX_inv  = inverse_spd(XX_new);
      matrix[50, 50] P_grid  = X_new * XX_inv * X_new';
      matrix[50, 50] I_G     = diag_matrix(rep_vector(1.0, 50));
      matrix[50, 50] IminP_g = I_G - P_grid;

      // Cross-cov (grid × train) for ARD GP with other coords fixed at 0
      matrix[50, N] K_cross;
      {
        real a2 = square(alpha);
        for (g in 1:50) {
          for (n in 1:N) {
            real r2 = 0;
            // vary j on grid; other predictors fixed at 0 (centered)
            r2 += square( (x_c_grid[g] - X_c[n, j]) / lambda[j] );
            for (m in 1:P) if (m != j)
              r2 += square( (0 - X_c[n, m]) / lambda[m] );
            K_cross[g, n] = a2 * exp(-0.5 * r2);
          }
        }
      }

      // Orthogonalize both sides consistently with model
      matrix[50, N] Ktilde_cross = IminP_g * K_cross * IminP_lin;

      // Partial dependence GP mean and full grid curve (ORIGINAL scale)
      f_mean_grid[, j]     = Ktilde_cross * a;
      linear_est_grid[, j] = rep_vector(beta0_orig, 50) + (sigma * theta[j]) * x_c_grid;
      full_est_grid[, j]   = linear_est_grid[, j] + f_mean_grid[, j];
    }
  }
}
