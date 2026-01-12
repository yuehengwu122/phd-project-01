data {
  int<lower=1> N;
  array[N] real x;                     // 1D inputs (raw)
  vector[N] y;                         // response (raw)
}

transformed data {
  // ----- Center x and y -----
  real x_bar = mean(x);
  real y_bar = mean(y);

  array[N] real x_c;
  vector[N] y_c;
  for (n in 1:N) x_c[n] = x[n] - x_bar;
  y_c = y - y_bar;

  // ----- Design for orthogonalization: X = [1, x_c] -----
  matrix[N, 2] X;
  for (n in 1:N) {
    X[n, 1] = 1.0;
    X[n, 2] = x_c[n];
  }
}

parameters {
  // --- GP hyperparameters (delta on sigma-scale) ---
  real<lower=0> lambda;        // length-scale
  real<lower=0> delta2;         // standardized GP amplitude

  // g-prior scale for slope (marginalized β1)
  real<lower=0> tau2;          // slope-only g-prior scale^2 (s=1)

  // noise
  real<lower=0> sigma;         // noise SD

  // intercept (centered scale)
  real mu;                     // mean of y_c
}

transformed parameters {
  real<lower=0> alpha = sqrt(delta2) * sigma;

  // Prior buffer (for bridge sampling)
  real lprior = 0;

  // half-t priors via truncated Student-t
  lprior += student_t_lpdf(delta2 | 1.15, 0, 7) - student_t_lccdf(0 | 1.15, 0, 7);
  lprior += inv_gamma_lpdf(lambda | 3, 0.8);

  lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);

  // mu ~ Student-t(3, 0, 10)
  lprior += student_t_lpdf(mu | 3, 0, 10);

  // g-prior on slope: tau2 ~ IG(1/2, 1/2 * s^2) with s = 1
  lprior += inv_gamma_lpdf(tau2 | 0.5, 0.5 * 1.0);
}

model {
  target += lprior;

  // ----- Projections for orthogonalization -----
  matrix[N,N] I = diag_matrix(rep_vector(1.0, N));
  matrix[2,2] XX = crossprod(X);
  matrix[2,N] S  = mdivide_left_spd(XX, X');
  matrix[N,N] P  = X * S;                         // projector onto col([1, x_c])
  matrix[N,N] I_minus_P = I - P;

  // slope-only projector P_s = x x' / (x' x), x = x_c
  vector[N] x_slope = col(X, 2);
  real xTx = dot_self(x_slope);
  row_vector[N] x_slope_row = to_row_vector(x_slope);
  matrix[N,N] P_s = (x_slope * x_slope_row) / xTx;

  // ----- Latent GP covariance (centered x) -----
  matrix[N, N] K = gp_exp_quad_cov(x_c, alpha, lambda);
  for (n in 1:N) K[n, n] += 1e-12;

  // ----- Orthogonalized covariance: Ktilde = (I - P) K (I - P), symmetrized -----
  matrix[N,N] Ktilde = 0.5 * (I_minus_P * K * I_minus_P
                              + (I_minus_P * K * I_minus_P)');

  // ----- Σ_y (CENTERED; ORTHOGONAL GP; slope marginalized) -----
  matrix[N, N] Sigma = Ktilde + square(sigma) * tau2 * P_s;
  for (n in 1:N) Sigma[n, n] += square(sigma);

  // ----- Mean on centered scale: intercept only (β1 marginalized) -----
  vector[N] mu_vec = rep_vector(mu, N);

  // Likelihood
  matrix[N, N] L = cholesky_decompose(Sigma);
  target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
}

generated quantities {
  // training-point estimates (ORIGINAL scale)
  vector[N] f_mean;         // E[f | y, θ] (orthogonal GP contribution)
  vector[N] full_est;       // (β0 + β1_mean x_c) + f_mean
  vector[N] linear_est;     // uses β1_mean

  // posterior summaries
  vector[2] beta_mean;      // [beta0_orig, beta1_mean]
  vector[2] theta_draw;     // θ = β / σ with β1_mean

  // 50-point grid (ORIGINAL scale)
  vector[50] x_grid;
  vector[50] f_mean_grid;
  vector[50] linear_est_grid;
  vector[50] full_est_grid;

  {
    // Projections (same as in model block)
    matrix[N,N] I = diag_matrix(rep_vector(1.0, N));
    matrix[2,2] XX = crossprod(X);
    matrix[2,N] S  = mdivide_left_spd(XX, X');
    matrix[N,N] P  = X * S;
    matrix[N,N] I_minus_P = I - P;

    vector[N] x_slope = col(X, 2);
    real xTx = dot_self(x_slope);
    row_vector[N] x_slope_row = to_row_vector(x_slope);
    matrix[N,N] P_s = (x_slope * x_slope_row) / xTx;

    // Latent GP covariance and orthogonalization
    matrix[N,N] K = gp_exp_quad_cov(x_c, alpha, lambda);
    for (n in 1:N) K[n,n] += 1e-12;

    matrix[N,N] Ktilde = 0.5 * (I_minus_P * K * I_minus_P
                                + (I_minus_P * K * I_minus_P)');

    // Σ_y and its Cholesky (centered scale)
    matrix[N,N] Sigma = Ktilde + square(sigma) * tau2 * P_s;
    for (n in 1:N) Sigma[n,n] += square(sigma);
    matrix[N,N] L = cholesky_decompose(Sigma);

    // a = Σ_y^{-1} (y_c - mu)
    vector[N] resid = y_c - rep_vector(mu, N);
    vector[N] v = mdivide_left_tri_low(L, resid);
    vector[N] a = mdivide_right_tri_low(v', L)';

    // Posterior mean of orthogonal GP at training inputs
    f_mean = Ktilde * a;        // Cov(f,y) Σ^{-1} resid  (with Ktilde)

    // Ridge (g-prior) posterior mean of slope β1 under Σ:
    real beta1_mean;
    {
      real t_shrink = tau2 / (1 + tau2);
      vector[N] y_tilde = resid - f_mean;   // remove GP & intercept
      real xTy = dot_product(x_slope, y_tilde);
      beta1_mean = t_shrink * (xTy / xTx);
    }

    // Intercept on ORIGINAL scale
    real beta0_orig = mu + y_bar;

    // Collect β means (original scale)
    beta_mean[1] = beta0_orig;
    beta_mean[2] = beta1_mean;

    // Linear trend on ORIGINAL scale
    linear_est = rep_vector(beta0_orig, N) + beta1_mean * col(X, 2);

    // θ = β / σ  (bookkeeping)
    theta_draw[1] = beta0_orig / sigma;
    theta_draw[2] = beta1_mean / sigma;

    // Total estimate on ORIGINAL scale
    full_est = linear_est + f_mean;

    // ---- Grid predictions (ORIGINAL scale) ----
    real xmin = min(x);
    real xmax = max(x);
    for (g in 1:50) x_grid[g] = xmin + (xmax - xmin) * (g - 1) / 49;

    // Center grid with same x mean
    real xbar = mean(x);
    matrix[50,2] X_new;
    for (g in 1:50) {
      real xc = x_grid[g] - xbar;
      X_new[g,1] = 1.0;
      X_new[g,2] = xc;
    }

    // Grid projector for orthogonalization
    matrix[50,50] I_G     = diag_matrix(rep_vector(1.0, 50));
    matrix[2,2]   XX_new  = crossprod(X_new);
    matrix[2,2]   XX_inv  = inverse_spd(XX_new);
    matrix[50,50] P_new   = X_new * XX_inv * X_new';
    matrix[50,50] IminP_n = I_G - P_new;

    // Cross-covariances (latent), then orthogonalize both sides
    matrix[50,N] K_cross;
    {
      real a2 = square(alpha);
      for (g in 1:50)
        for (n in 1:N) {
          real xc_g = X_new[g,2];
          real xc_n = X[n,2];
          real r = (xc_g - xc_n) / lambda;
          K_cross[g,n] = a2 * exp(-0.5 * r * r);
        }
    }
    matrix[50,N] Ktilde_cross = IminP_n * K_cross * (I - P);

    // Grid means (ORIGINAL scale)
    f_mean_grid      = Ktilde_cross * a;
    linear_est_grid  = rep_vector(beta0_orig, 50) + beta1_mean * X_new[,2];
    full_est_grid    = linear_est_grid + f_mean_grid;
  }
}
