data {
  int<lower=1> N;
  array[N] real x;                     // 1D inputs (raw)
  vector[N] y;                         // response (raw)
  array[N] int<lower=0,upper=1> d;     // group indicator (0/1)
}

transformed data {
  // ----- Center x and y -----
  real x_bar = mean(x);
  real y_bar = mean(y);

  array[N] real x_c;
  vector[N]    y_c;
  for (n in 1:N) x_c[n] = x[n] - x_bar;
  y_c = y - y_bar;

  // ----- Design for orthogonalization: X = [1, x_c] -----
  matrix[N, 2] X;
  for (n in 1:N) {
    X[n, 1] = 1.0;
    X[n, 2] = x_c[n];
  }

  // ----- Centered group dummy: z = d - 0.5 -----
  vector[N] z;
  for (n in 1:N) z[n] = d[n] - 0.5;
}

parameters {
  // --- GP hyperparameters (delta on sigma-scale) ---
  // common GP g(x)
  real<lower=0> lambda_g;      // length-scale for g
  real<lower=0> delta_g;       // standardized GP amplitude for g

  // interaction GP h(x)
  real<lower=0> lambda_h;      // length-scale for h
  real<lower=0> delta_h;       // standardized GP amplitude for h

  // g-prior scale for slope (marginalized β1)
  real<lower=0> tau2;          // slope-only g-prior scale^2 (s=1)

  // noise
  real<lower=0> sigma;         // noise SD

  // intercept (centered scale)
  real mu;                     // mean of y_c
}

transformed parameters {
  // Convert standardized amplitudes to marginal SDs
  real<lower=0> alpha_g = delta_g * sigma;
  real<lower=0> alpha_h = delta_h * sigma;

  // Prior buffer (for bridge sampling)
  real lprior = 0;

  // half-t priors via truncated Student-t
  lprior += student_t_lpdf(delta_g | 4, 0, 2.7) - student_t_lccdf(0 | 4, 0, 2.7);
  lprior += inv_gamma_lpdf(lambda_g | 3, 0.8);

  lprior += student_t_lpdf(delta_h | 4, 0, 2.7) - student_t_lccdf(0 | 4, 0, 2.7);
  lprior += inv_gamma_lpdf(lambda_h | 3, 0.8);

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

  // ----- Latent GP covariances (centered x) -----
  matrix[N, N] Kg = gp_exp_quad_cov(x_c, alpha_g, lambda_g);
  matrix[N, N] Kh = gp_exp_quad_cov(x_c, alpha_h, lambda_h);
  for (n in 1:N) { Kg[n, n] += 1e-12; Kh[n, n] += 1e-12; }

  // ----- Orthogonalized covariances: Ktilde = (I - P) K (I - P) -----
  matrix[N,N] Ktilde_g = 0.5 * (I_minus_P * Kg * I_minus_P
                                + (I_minus_P * Kg * I_minus_P)');
  matrix[N,N] Ktilde_h = 0.5 * (I_minus_P * Kh * I_minus_P
                                + (I_minus_P * Kh * I_minus_P)');

  // ----- Centered group-diagonal -----
  matrix[N,N] Dz = diag_matrix(z);

  // ----- Σ_y (CENTERED scale; ORTHOGONAL GPs; slope marginalized) -----
  matrix[N, N] Sigma = Ktilde_g + Dz * Ktilde_h * Dz + square(sigma) * tau2 * P_s;
  for (n in 1:N) Sigma[n, n] += square(sigma);

  // ----- Mean on centered scale: intercept only (β1 marginalized) -----
  vector[N] mu_vec = rep_vector(mu, N);

  // Likelihood
  matrix[N, N] L = cholesky_decompose(Sigma);
  target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
}

generated quantities {
  // training-point trends (ORIGINAL scale)
  vector[N] linear_est;        // (mu + y_bar) + beta1_mean * x_c  (β1 from ridge mean)
  vector[N] g_mean;            // E[g | y, θ]  (orthogonal GP)
  vector[N] h_mean;            // E[h | y, θ]  (orthogonal GP)
  vector[N] trend_common;      // lin + g
  vector[N] trend_group0;      // lin + g - 0.5 h
  vector[N] trend_group1;      // lin + g + 0.5 h

  // posterior summaries
  vector[2] beta_mean;         // [beta0_orig, beta1_mean]
  vector[2] theta_draw;        // θ = β / σ  (with β1_mean)

  // 50-point grid (ORIGINAL scale)
  vector[50] x_grid;
  vector[50] linear_est_grid;
  vector[50] g_mean_grid;
  vector[50] h_mean_grid;
  vector[50] trend_common_grid;
  vector[50] trend_group0_grid;
  vector[50] trend_group1_grid;

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

    // Latent GP covariances and orthogonalization
    matrix[N,N] Kg = gp_exp_quad_cov(x_c, delta_g * sigma, lambda_g);
    matrix[N,N] Kh = gp_exp_quad_cov(x_c, delta_h * sigma, lambda_h);
    for (n in 1:N) { Kg[n,n] += 1e-12; Kh[n,n] += 1e-12; }

    matrix[N,N] Ktilde_g = 0.5 * (I_minus_P * Kg * I_minus_P
                                  + (I_minus_P * Kg * I_minus_P)');
    matrix[N,N] Ktilde_h = 0.5 * (I_minus_P * Kh * I_minus_P
                                  + (I_minus_P * Kh * I_minus_P)');

    // Σ_y and its Cholesky (centered scale)
    matrix[N,N] Dz = diag_matrix(z);
    matrix[N,N] Sigma = Ktilde_g + Dz * Ktilde_h * Dz + square(sigma) * tau2 * P_s;
    for (n in 1:N) Sigma[n,n] += square(sigma);
    matrix[N,N] L = cholesky_decompose(Sigma);

    // a = Σ_y^{-1} (y_c - mu)
    vector[N] resid = y_c - rep_vector(mu, N);
    vector[N] v = mdivide_left_tri_low(L, resid);
    vector[N] a = mdivide_right_tri_low(v', L)';

    // Posterior means of orthogonal GPs at training inputs
    g_mean = Ktilde_g * a;           // Cov(g,y) Σ^{-1} resid
    h_mean = Ktilde_h * (Dz * a);    // Cov(h,y) = Ktilde_h Dz

    // Ridge (g-prior) posterior mean of slope β1 under Σ:
    // shrinkage t = tau2 / (1 + tau2), using y_tilde with g removed
    real beta1_mean;
    {
      real t_shrink = tau2 / (1 + tau2);
      vector[N] y_tilde = resid - g_mean;   // remove common GP & intercept
      real xTy = dot_product(x_slope, y_tilde);
      beta1_mean = t_shrink * (xTy / xTx);
    }

    // Intercept on ORIGINAL scale: add back y_bar
    real beta0_orig = mu + y_bar;

    // Collect β means (original scale)
    beta_mean[1] = beta0_orig;
    beta_mean[2] = beta1_mean;

    // Linear trend on ORIGINAL scale
    linear_est = rep_vector(beta0_orig, N) + beta1_mean * col(X, 2);

    // θ = β / σ  (for SDR/KDE bookkeeping)
    theta_draw[1] = beta0_orig / sigma;
    theta_draw[2] = beta1_mean / sigma;

    // Total trends on ORIGINAL scale
    trend_common = linear_est + g_mean;
    trend_group0 = trend_common - 0.5 * h_mean;
    trend_group1 = trend_common + 0.5 * h_mean;

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
    matrix[50,N] K_no_g;
    matrix[50,N] K_no_h;
    {
      real a2_g = square(delta_g * sigma);
      real a2_h = square(delta_h * sigma);
      for (g in 1:50)
        for (n in 1:N) {
          real xc_g = X_new[g,2];
          real xc_n = X[n,2];
          real r_g = (xc_g - xc_n) / lambda_g;
          real r_h = (xc_g - xc_n) / lambda_h;
          K_no_g[g,n] = a2_g * exp(-0.5 * r_g * r_g);
          K_no_h[g,n] = a2_h * exp(-0.5 * r_h * r_h);
        }
    }
    matrix[50,N] Ktilde_no_g = IminP_n * K_no_g * (I - P);
    matrix[50,N] Ktilde_no_h = IminP_n * K_no_h * (I - P);

    // Grid means (ORIGINAL scale)
    vector[50] lin_grid = rep_vector(beta0_orig, 50) + beta1_mean * X_new[,2];
    g_mean_grid        = Ktilde_no_g * a;
    h_mean_grid        = Ktilde_no_h * (Dz * a);
    linear_est_grid    = lin_grid;
    trend_common_grid  = lin_grid + g_mean_grid;
    trend_group0_grid  = trend_common_grid - 0.5 * h_mean_grid;
    trend_group1_grid  = trend_common_grid + 0.5 * h_mean_grid;
  }
}
