functions {
  /**
   * Unit-amplitude SE kernel from precomputed squared distances.
   *
   * R[i,j] = exp(-0.5 * D2[i,j] / l^2)  for i != j
   * R[i,i] = 1 + nugget
   */
  matrix se_kernel_unit(matrix D2, real l, real nugget) {
    int N = rows(D2);
    real inv_l2 = inv(square(l));
    matrix[N, N] R = exp(-0.5 * inv_l2 * D2);
    for (n in 1:N)
      R[n, n] = 1.0 + nugget;
    return R;
  }
}

data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;
  vector[N] y;

  int<lower=0> d_order;
  real<lower=0> kappa_delta;
  int<lower=1> df_delta;
  real<lower=0> s_delta;
  real logC_base;
}

transformed data {
  real ig_a = 2.0;
  real ig_b = 0.5 * N;

  real y_bar = mean(y);
  vector[N] y_c = y - y_bar;

  row_vector[P] x_bar;
  matrix[N, P] X_c;
  for (j in 1:P) {
    x_bar[j] = mean(col(X, j));
    X_c[, j] = col(X, j) - x_bar[j];
  }

  matrix[P, P] XtX = crossprod(X_c);
  matrix[P, P] XtX_ridge = XtX + 1e-10 * diag_matrix(rep_vector(1.0, P));
  matrix[P, P] V_theta_base = inverse_spd(XtX_ridge);
  matrix[P, P] L_V_theta_base = cholesky_decompose(V_theta_base);

  int R_cols = P + 1;
  matrix[N, R_cols] Z_aug;
  Z_aug[, 1] = rep_vector(1.0, N);
  Z_aug[, 2:R_cols] = X_c;
  matrix[N, R_cols] Q_Z = qr_thin_Q(Z_aug);

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
  vector[P] beta1 = sigma * theta;

  real lprior = 0;

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

  for (j in 1:P)
    lprior += inv_gamma_lpdf(lambda[j] | 3, 0.8);

  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
            - student_t_lccdf(0 | 3, 0, 2.5);

  {
    real g_theta = exp(log_g_theta);
    lprior += inv_gamma_lpdf(g_theta | ig_a, ig_b) + log_g_theta;
  }
}

model {
  target += lprior;
  z_theta ~ normal(0, 1);

  int R = R_cols;
  matrix[N, N] S_R = rep_matrix(0.0, N, N);
  matrix[N, R] S_U = rep_matrix(0.0, N, R);
  matrix[R, R] S_M = rep_matrix(0.0, R, R);

  for (j in 1:P) {
    matrix[N, N] Rj = se_kernel_unit(D2[j], lambda[j], 1e-8);
    matrix[N, R] Uj = Rj * Q_Z;
    matrix[R, R] Mj = Q_Z' * Uj;

    // trace(P_perp * Rj * P_perp) = trace(Rj) - trace(Q_Z' * Rj * Q_Z)
    real tr_proj = trace(Rj) - trace(Mj);
    real c_j = tr_proj / N;
    real c_safe = fmax(c_j, 1e-10);
    real w_j = square(sigma * delta[j]) / c_safe;

    S_R += w_j * Rj;
    S_U += w_j * Uj;
    S_M += w_j * Mj;
  }

  matrix[N, N] Sigma = S_R
                        - S_U * Q_Z'
                        - Q_Z * S_U'
                        + Q_Z * S_M * Q_Z';

  {
    real s2 = square(sigma);
    for (n in 1:N)
      Sigma[n, n] += s2;
  }

  Sigma = 0.5 * (Sigma + Sigma');

  vector[N] mu_vec = X_c * beta1;
  {
    matrix[N, N] L = cholesky_decompose(Sigma);
    target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
  }
}
