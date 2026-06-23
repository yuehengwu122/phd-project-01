// =============================================================================
// fit_interaction_only_model.stan
//
// Reduced model for bridge-sampling validation of the main nonlinear BF.
// This is the full interaction model with delta_main = 0 (f_main removed):
//
//   y_i = beta_0 + beta_1 x_i + z_i f_int(x_i) + eps_i
//   f_int ~ GP(0, sigma^2 delta_int^2 K~^*_{lambda_int})
//   eps_i ~ N(0, sigma^2)
//
// K~^*_int is orthogonalized against span{1, x_c} and trace-normalized so
// delta_int measures the projected interaction SD on the response scale.
// =============================================================================

functions {
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
  vector[N] x;
  vector[N] y;
  vector[N] z;

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
  real x_bar = mean(x);
  vector[N] y_c = y - y_bar;
  vector[N] x_c = x - x_bar;

  real xtx = dot_self(x_c) + 1e-10;
  real V_theta_base = inv(xtx);
  real sd_theta_base = sqrt(V_theta_base);

  int R_cols = 2;
  matrix[N, R_cols] Z_aug;
  Z_aug[, 1] = rep_vector(1.0, N);
  Z_aug[, 2] = x_c;
  matrix[N, R_cols] Q_Z = qr_thin_Q(Z_aug);

  matrix[N, N] D2 = rep_matrix(0.0, N, N);
  if (N > 1) {
    for (n in 1:(N - 1)) {
      for (m in (n + 1):N) {
        real d = x_c[n] - x_c[m];
        real d2 = d * d;
        D2[n, m] = d2;
        D2[m, n] = d2;
      }
    }
  }
}

parameters {
  real log_g_theta;
  real z_theta;

  real<lower=0> lambda_int;
  real<lower=1e-12> delta_int;

  real<lower=0> sigma;
}

transformed parameters {
  real theta = exp(0.5 * log_g_theta) * sd_theta_base * z_theta;
  real beta1 = sigma * theta;

  real lprior = 0;

  {
    real s_eff = s_delta * kappa_delta;
    lprior += student_t_lpdf(delta_int | df_delta, 0, s_eff)
              - student_t_lccdf(0 | df_delta, 0, s_eff);
    if (d_order > 0)
      lprior += d_order * log(delta_int);
    lprior -= (logC_base + d_order * log(kappa_delta));
  }

  lprior += inv_gamma_lpdf(lambda_int | 3, 0.8);

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

  matrix[N, N] R_int = se_kernel_unit(D2, lambda_int, 1e-8);
  matrix[N, R_cols] U_int = R_int * Q_Z;
  matrix[R_cols, R_cols] M_int = Q_Z' * U_int;
  matrix[N, N] Ktilde_int = R_int
                            - U_int * Q_Z'
                            - Q_Z * U_int'
                            + Q_Z * M_int * Q_Z';

  real tr_int = 0;
  for (n in 1:N)
    tr_int += square(z[n]) * Ktilde_int[n, n];
  real c_int = tr_int / N;
  real c_int_safe = fmax(c_int, 1e-10);
  real w_int = square(sigma * delta_int) / c_int_safe;

  matrix[N, N] Sigma;
  for (i in 1:N)
    for (j in 1:N)
      Sigma[i, j] = w_int * z[i] * Ktilde_int[i, j] * z[j];

  {
    real s2 = square(sigma);
    for (n in 1:N)
      Sigma[n, n] += s2;
  }

  Sigma = 0.5 * (Sigma + Sigma');

  vector[N] mu_vec = x_c * beta1;
  {
    matrix[N, N] L = cholesky_decompose(Sigma);
    target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
  }
}
