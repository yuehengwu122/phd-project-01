functions {
  /**
   * Build SE kernel from precomputed squared-distance matrix D2.
   * Uses vectorized operations to minimize autodiff tape nodes.
   *
   * K[i,j] = alpha^2 * exp(-0.5 * D2[i,j] / l^2)  for i != j
   * K[i,i] = alpha^2 + nugget
   */
  matrix se_kernel(matrix D2, real alpha, real l, real nugget) {
    int N = rows(D2);
    real a2 = square(alpha);
    real inv_l2 = inv(square(l));
    // Vectorized: exp operates element-wise on the whole matrix at once
    // This produces ONE autodiff node for the entire exp(), not N^2 nodes.
    matrix[N, N] K = a2 * exp(-0.5 * inv_l2 * D2);
    // Fix diagonal (D2 diagonal is 0, so exp(0)=1, giving a2; add nugget)
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

  // ----- Full-rank orthogonal projector against [1 | X_c] -----
  //
  // Design matrix augmented with intercept column:
  //   Z = [1_N , X_c]   (N x (P+1))
  //
  // Orthogonal complement projector:
  //   P_perp = I_N - Z (Z'Z)^{-1} Z'
  //
  // We precompute the QR of Z so that P_perp can be applied as:
  //   P_perp * A = A - Q_Z * (Q_Z' * A)
  //
  // This avoids ever forming the dense N x N projector.
  // Q_Z is N x (P+1), which is much smaller than N x N.

  int R_cols = P + 1;  // rank of [1 | X_c]
  matrix[N, R_cols] Z_aug;
  Z_aug[, 1] = rep_vector(1.0, N);
  Z_aug[, 2:R_cols] = X_c;

  // Thin QR: Z_aug = Q_Z * R_Z, Q_Z is N x (P+1) orthonormal
  matrix[N, R_cols] Q_Z = qr_thin_Q(Z_aug);

  // ----- Precompute squared distances for each predictor -----
  array[P] matrix[N, N] D2;
  for (j in 1:P) {
    D2[j] = rep_matrix(0.0, N, N);
    if (N > 1) {
      for (n in 1:(N - 1)) {
        for (m in (n + 1):N) {
          real d = X_c[n, j] - X_c[m, j];
          real d2 = d * d;
          D2[j][n, m] = d2;
          D2[j][m, n] = d2;
        }
      }
    }
  }
}

parameters {
  real log_g_theta;
  vector[P] z_theta;
  vector<lower=0>[P] lambda;
  vector<lower=1e-12>[P] delta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[P] theta = exp(0.5 * log_g_theta) * (L_V_theta_base * z_theta);
  vector<lower=0>[P] alpha = delta * sigma;
  vector[P] beta1 = sigma * theta;

  real lprior = 0;

  // ---- moment prior on delta ----
  {
    real s_eff = s_delta * kappa_delta;
    for (j in 1:P) {
      real delt = delta[j];
      lprior += student_t_lpdf(delt | df_delta, 0, s_eff)
                - student_t_lccdf(0 | df_delta, 0, s_eff);
      if (d_order > 0)
        lprior += d_order * log(delt);
      lprior -= (logC_base + d_order * log(kappa_delta));
    }
  }

  // lambda prior
  for (j in 1:P)
    lprior += inv_gamma_lpdf(lambda[j] | 3, 0.8);

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

  // ----- Build Sigma = sum_j Ktilde_j + sigma^2 * I -----
  //
  // For each predictor j:
  //   K_j         = SE kernel (depends on alpha[j], lambda[j])
  //   Ktilde_j    = P_perp * K_j * P_perp
  //
  // where P_perp = I - Q_Z * Q_Z'.
  //
  // Expanding:
  //   Ktilde_j = K_j - Q_Z * (Q_Z' * K_j) - (K_j * Q_Z) * Q_Z' + Q_Z * (Q_Z' * K_j * Q_Z) * Q_Z'
  //
  // Key optimization: rather than forming four N x N temporaries,
  // we accumulate everything into Sigma as follows:
  //
  //   Let U_j = K_j * Q_Z              (N x R)  -- one N x N times N x R multiply
  //   Let M_j = Q_Z' * U_j             (R x R)  -- one R x N times N x R multiply
  //
  //   Then:  Ktilde_j = K_j - U_j * Q_Z' - Q_Z * U_j' + Q_Z * M_j * Q_Z'
  //
  // We accumulate across all j:
  //   Sigma = sum_j K_j  -  (sum_j U_j) * Q_Z'  -  Q_Z * (sum_j U_j)'  +  Q_Z * (sum_j M_j) * Q_Z'  +  sigma^2 * I
  //
  // This means we need:
  //   S_K   = sum_j K_j          (N x N)  -- unavoidable, one N x N per j
  //   S_U   = sum_j U_j          (N x R)  -- accumulated cheaply
  //   S_M   = sum_j M_j          (R x R)  -- accumulated cheaply
  //
  // Then ONE final assembly:
  //   Sigma = S_K - S_U * Q_Z' - Q_Z * S_U' + Q_Z * S_M * Q_Z' + sigma^2 * I
  //
  // This replaces P * 7 dense N x N operations with:
  //   P kernel constructions + P cheap (N x R) multiplies + ONE final N x N assembly.

  int R = R_cols;  // = P + 1

  matrix[N, N] S_K = rep_matrix(0.0, N, N);
  matrix[N, R] S_U = rep_matrix(0.0, N, R);
  matrix[R, R] S_M = rep_matrix(0.0, R, R);

  for (j in 1:P) {
    matrix[N, N] Kj = se_kernel(D2[j], alpha[j], lambda[j], 1e-8);
    matrix[N, R] Uj = Kj * Q_Z;       // N x R  (the only big multiply per j)
    matrix[R, R] Mj = Q_Z' * Uj;      // R x R

    S_K += Kj;
    S_U += Uj;
    S_M += Mj;
  }

  // Final assembly: Sigma = S_K - S_U * Q_Z' - Q_Z * S_U' + Q_Z * S_M * Q_Z' + sigma^2 I
  matrix[N, N] Sigma = S_K
                        - S_U * Q_Z'
                        - Q_Z * S_U'
                        + Q_Z * S_M * Q_Z';

  // Add noise variance
  {
    real s2 = square(sigma);
    for (n in 1:N)
      Sigma[n, n] += s2;
  }

  // Symmetrize (guard against floating-point asymmetry)
  Sigma = 0.5 * (Sigma + Sigma');

  // Likelihood
  vector[N] mu_vec = X_c * beta1;
  {
    matrix[N, N] L = cholesky_decompose(Sigma);
    target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
  }
}
