functions {
  /**
   * Build SE kernel from precomputed squared-distance matrix D2.
   */
  matrix se_kernel(matrix D2, real alpha, real l, real nugget) {
    int N = rows(D2);
    real a2 = square(alpha);
    real inv_l2 = inv(square(l));
    matrix[N, N] K = a2 * exp(-0.5 * inv_l2 * D2);
    for (n in 1:N)
      K[n, n] = a2 + nugget;
    return K;
  }
}

data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;
  vector[N] y;

  // ---- Selection masks (0/1). Length P for both. ----
  // include_linear[j]    = 1  => predictor j enters the linear part
  // include_nonlinear[j] = 1  => predictor j enters the additive GP sum
  array[P] int<lower=0, upper=1> include_linear;
  array[P] int<lower=0, upper=1> include_nonlinear;

  // ---- moment-prior controls (computed outside Stan) ----
  int<lower=0> d_order;
  real<lower=0> kappa_delta;
  int<lower=1> df_delta;
  real<lower=0> s_delta;
  real logC_base;
}

transformed data {
  real ig_a = 2.0;
  real ig_b = 0.5 * N;

  // ----- Counts of active predictors -----
  int P_lin = 0;
  int P_nl  = 0;
  for (j in 1:P) {
    P_lin += include_linear[j];
    P_nl  += include_nonlinear[j];
  }
  // We require at least one linear predictor so that the projector
  // is well-defined (the intercept column alone is also fine, but
  // keeping P_lin >= 0 via the `if` branches below handles P_lin = 0).
  // P_nl = 0 is allowed (model reduces to linear-only), but for the
  // intended use case (leave-one-nonlinear-out) we expect P_nl >= 1.

  // ----- Index sets for active predictors -----
  array[P_lin] int idx_lin;
  array[P_nl]  int idx_nl;
  {
    int a = 1;
    int b = 1;
    for (j in 1:P) {
      if (include_linear[j] == 1)    { idx_lin[a] = j; a += 1; }
      if (include_nonlinear[j] == 1) { idx_nl[b]  = j; b += 1; }
    }
  }

  // ----- Center y and each column of X -----
  real y_bar = mean(y);
  vector[N] y_c = y - y_bar;

  row_vector[P] x_bar;
  matrix[N, P] X_c;
  for (j in 1:P) {
    x_bar[j] = mean(col(X, j));
    X_c[, j] = col(X, j) - x_bar[j];
  }

  // ----- Reduced linear design matrix (only linear-included columns) -----
  matrix[N, P_lin] X_lin;
  for (a in 1:P_lin)
    X_lin[, a] = col(X_c, idx_lin[a]);

  // ---- mixture-of-g prior base covariance, on the REDUCED linear block ----
  // If P_lin == 0 we still need a placeholder of size 0x0 which Stan handles.
  matrix[P_lin, P_lin] XtX_lin = crossprod(X_lin);
  matrix[P_lin, P_lin] XtX_lin_ridge =
      XtX_lin + 1e-10 * diag_matrix(rep_vector(1.0, P_lin));
  matrix[P_lin, P_lin] V_theta_base = inverse_spd(XtX_lin_ridge);
  matrix[P_lin, P_lin] L_V_theta_base = cholesky_decompose(V_theta_base);

  // ----- Orthogonal projector against [1 | X_c] (FULL X_c) -----
  //
  // IMPORTANT: the projector is built from ALL P linear columns,
  // independently of include_linear. This fixes the function space
  // of each GP component f_j to be the orthogonal complement of the
  // FULL linear span. Without this invariance, delta_j would mean
  // different things in different sub-models (nonlinear-only in some,
  // linear+nonlinear in others) and the Bayes factor for the nonlinear
  // component would be ill-defined. The projector defines what
  // "nonlinear" means, and that meaning must be invariant across the
  // models being compared.
  int R_cols = P + 1;
  matrix[N, R_cols] Z_aug;
  Z_aug[, 1] = rep_vector(1.0, N);
  Z_aug[, 2:R_cols] = X_c;

  matrix[N, R_cols] Q_Z = qr_thin_Q(Z_aug);

  // ----- Precompute squared distances ONLY for nonlinear-active predictors -----
  array[P_nl] matrix[N, N] D2;
  for (b in 1:P_nl) {
    int j = idx_nl[b];
    D2[b] = rep_matrix(0.0, N, N);
    if (N > 1) {
      for (n in 1:(N - 1)) {
        for (m in (n + 1):N) {
          real d = X_c[n, j] - X_c[m, j];
          real d2 = d * d;
          D2[b][n, m] = d2;
          D2[b][m, n] = d2;
        }
      }
    }
  }
}

parameters {
  real log_g_theta;
  vector[P_lin] z_theta;
  vector<lower=0>[P_nl] lambda;
  vector<lower=1e-12>[P_nl] delta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[P_lin] theta = exp(0.5 * log_g_theta) * (L_V_theta_base * z_theta);
  vector<lower=0>[P_nl] alpha = delta * sigma;
  vector[P_lin] beta1_lin = sigma * theta;

  real lprior = 0;

  // ---- moment prior on delta (only over active nonlinear components) ----
  {
    real s_eff = s_delta * kappa_delta;
    for (b in 1:P_nl) {
      real delt = delta[b];
      lprior += student_t_lpdf(delt | df_delta, 0, s_eff)
                - student_t_lccdf(0 | df_delta, 0, s_eff);
      if (d_order > 0)
        lprior += d_order * log(delt);
      lprior -= (logC_base + d_order * log(kappa_delta));
    }
  }

  // lambda prior (only over active nonlinear components)
  for (b in 1:P_nl)
    lprior += inv_gamma_lpdf(lambda[b] | 3, 0.8);

  // sigma prior
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
            - student_t_lccdf(0 | 3, 0, 2.5);

  // ---- IG prior on g_theta = exp(log_g_theta) ----
  {
    real g_theta = exp(log_g_theta);
    lprior += inv_gamma_lpdf(g_theta | ig_a, ig_b) + log_g_theta;
  }
}

model {
  target += lprior;
  z_theta ~ normal(0, 1);

  // ----- Build Sigma = sum_{b in active} Ktilde_b + sigma^2 * I -----
  int R = R_cols;

  matrix[N, N] S_K = rep_matrix(0.0, N, N);
  matrix[N, R] S_U = rep_matrix(0.0, N, R);
  matrix[R, R] S_M = rep_matrix(0.0, R, R);

  for (b in 1:P_nl) {
    matrix[N, N] Kj = se_kernel(D2[b], alpha[b], lambda[b], 1e-8);
    matrix[N, R] Uj = Kj * Q_Z;
    matrix[R, R] Mj = Q_Z' * Uj;

    S_K += Kj;
    S_U += Uj;
    S_M += Mj;
  }

  matrix[N, N] Sigma = S_K
                       - S_U * Q_Z'
                       - Q_Z * S_U'
                       + Q_Z * S_M * Q_Z';

  {
    real s2 = square(sigma);
    for (n in 1:N)
      Sigma[n, n] += s2;
  }

  Sigma = 0.5 * (Sigma + Sigma');

  // Mean: only the active linear columns contribute
  vector[N] mu_vec;
  if (P_lin > 0)
    mu_vec = X_lin * beta1_lin;
  else
    mu_vec = rep_vector(0.0, N);

  {
    matrix[N, N] L = cholesky_decompose(Sigma);
    target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
  }
}

generated quantities {
  // ---- Full-length beta1 and delta for reporting convenience ----
  // Excluded entries are set to 0 so downstream code can still index by
  // the original predictor position.
  vector[P] beta1 = rep_vector(0.0, P);
  vector[P] delta_full = rep_vector(0.0, P);
  vector[P] lambda_full = rep_vector(0.0, P);
  for (a in 1:P_lin)
    beta1[idx_lin[a]] = beta1_lin[a];
  for (b in 1:P_nl) {
    delta_full[idx_nl[b]]  = delta[b];
    lambda_full[idx_nl[b]] = lambda[b];
  }
}
