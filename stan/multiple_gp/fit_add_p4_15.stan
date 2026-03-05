functions {
  // Build SE kernel matrix from precomputed squared distances D2
  // K[n,m] = alpha^2 * exp(-0.5 * D2[n,m] / l^2)
  matrix se_cov_from_D2(matrix D2, real alpha, real l) {
    int N = rows(D2);
    matrix[N, N] K = rep_matrix(0.0, N, N);
    real a2 = square(alpha);
    real inv_l2 = inv(square(l));

    // diagonal
    for (n in 1:N)
      K[n, n] = a2;

    // upper triangle
    if (N > 1) {
      for (n in 1:(N - 1)) {
        for (m in (n + 1):N) {
          real k_nm = a2 * exp(-0.5 * D2[n, m] * inv_l2);
          K[n, m] = k_nm;
          K[m, n] = k_nm;
        }
      }
    }
    return K;
  }
}

data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;
  vector[N] y;
}

transformed data {
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

  // ----- Precompute rank-2 projector ingredients for each predictor -----
  array[P] matrix[N, 2] Q_tr;
  array[P] matrix[2, 2] R_tr;

  for (j in 1:P) {
    Q_tr[j][, 1] = rep_vector(1.0, N);
    Q_tr[j][, 2] = X_c[, j];
    {
      matrix[2, 2] G = crossprod(Q_tr[j]);   // Q'Q
      R_tr[j] = inverse_spd(G);
    }
  }

  // ----- Precompute squared distances for each predictor -----
  array[P] matrix[N, N] D2;
  for (j in 1:P) {
    // initialize to zeros (safe)
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
  // mixture-of-g prior scale for theta (log-scale)
  real log_g_theta;

  vector[P] z_theta; // standard normal

  vector<lower=0>[P] lambda;
  vector<lower=0>[P] delta4;

  real<lower=0> sigma;
}

transformed parameters {
  // theta | g_theta ~ N(0, g_theta * (X'X)^{-1})
  real<lower=0> g_theta = exp(log_g_theta);
  vector[P] theta = exp(0.5 * log_g_theta) * (L_V_theta_base * z_theta);

  // GP marginal SD per predictor
  vector<lower=0>[P] delta = delta4 ^ (1.0 / 4.0);
  vector<lower=0>[P] alpha = delta * sigma;

  // linear slopes on centered-x scale
  vector[P] beta1 = sigma * theta;

  real lprior = 0;

  for (j in 1:P)
    lprior += student_t_lpdf(delta4[j] | 0.7, 0, 25) - student_t_lccdf(0 | 0.7, 0, 25);

  for (j in 1:P)
    lprior += inv_gamma_lpdf(lambda[j] | 3, 0.8);

  lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);

  // log_g_theta ~ Normal(log(N), 0.5)
  lprior += normal_lpdf(log_g_theta | log(N), 0.5);
}

model {
  target += lprior;
  z_theta ~ normal(0, 1);

  // ----- Build Sigma = sum_j Ktilde_j + sigma^2 I -----
  matrix[N, N] Sigma = rep_matrix(0.0, N, N);

  for (j in 1:P) {
    // Kj from precomputed distances (N x N)
    matrix[N, N] Kj = se_cov_from_D2(D2[j], alpha[j], lambda[j]);
    for (n in 1:N) Kj[n, n] += 1e-12;

    // Ktilde = K - PK - KP + PKP, where P = Q R Q'
    matrix[N, 2] Q = Q_tr[j];
    matrix[2, 2] R = R_tr[j];

    // U = K Q  (N x 2)
    matrix[N, 2] U = Kj * Q;

    // PK = Q * (R * U')  (N x N)
    matrix[2, N] RU_t = R * U';
    matrix[N, N] PK = Q * RU_t;

    // M = Q' K Q = Q' U (2 x 2)
    matrix[2, 2] M = Q' * U;

    // PKP = Q * (R M R) * Q'
    matrix[2, 2] B = R * M * R;
    matrix[N, N] PKP = Q * B * Q';

    matrix[N, N] Ktilde = Kj - PK - PK' + PKP;
    Ktilde = 0.5 * (Ktilde + Ktilde'); // symmetrize

    Sigma += Ktilde;
  }

  for (n in 1:N)
    Sigma[n, n] += square(sigma);

  // Mean on centered scale: linear part
  vector[N] mu_vec = X_c * beta1;

  {
    matrix[N, N] L = cholesky_decompose(Sigma);
    target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
  }
}
