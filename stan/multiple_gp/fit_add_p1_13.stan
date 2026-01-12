functions {
  // Build SE kernel matrix from precomputed squared distances D2
  // K[n,m] = alpha^2 * exp(-0.5 * D2[n,m] / l^2)
  matrix se_cov_from_D2(matrix D2, real alpha, real l) {
    int N = rows(D2);
    matrix[N, N] K;
    real a2 = square(alpha);
    real inv_l2 = inv(square(l));  // 1 / l^2

    for (n in 1:N) {
      K[n, n] = a2;  // exp(0)=1
      for (m in (n + 1):N) {
        real k_nm = a2 * exp(-0.5 * D2[n, m] * inv_l2);
        K[n, m] = k_nm;
        K[m, n] = k_nm;
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

  // ----- Precompute rank-2 projector ingredients for each predictor -----
  array[P] matrix[N, 2] Q_tr;
  array[P] matrix[2, 2] R_tr;
  vector[P] xjxj;

  for (j in 1:P) {
    Q_tr[j][, 1] = rep_vector(1.0, N);
    Q_tr[j][, 2] = X_c[, j];

    {
      matrix[2, 2] G = crossprod(Q_tr[j]);   // Q'Q
      R_tr[j] = inverse_spd(G);
    }

    xjxj[j] = dot_self(X_c[, j]);
  }

  // ----- Precompute squared distances for each predictor -----
  array[P] matrix[N, N] D2;
  for (j in 1:P) {
    for (n in 1:N) {
      D2[j][n, n] = 0;
      for (m in (n + 1):N) {
        real d = X_c[n, j] - X_c[m, j];
        real d2 = d * d;
        D2[j][n, m] = d2;
        D2[j][m, n] = d2;
      }
    }
  }
}

parameters {
  vector<lower=0>[P] lambda;
  vector<lower=0>[P] delta;

  vector[P] theta;

  real<lower=0> sigma;
  real mu;
}

transformed parameters {
  vector<lower=0>[P] alpha = delta * sigma;
  vector[P] beta1 = sigma * theta;

  real lprior = 0;

  // delta_j ~ half-t(4, 2.7)
  for (j in 1:P)
    lprior += student_t_lpdf(delta[j] | 4, 0, 2.7) - student_t_lccdf(0 | 4, 0, 2.7);

  // lambda_j ~ IG(3, 0.8)
  for (j in 1:P)
    lprior += inv_gamma_lpdf(lambda[j] | 3, 0.8);

  // sigma ~ half-t(3, 2.5)
  lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);

  // mu ~ t(3, 0, 10)
  lprior += student_t_lpdf(mu | 3, 0, 10);

  // theta_j ~ Cauchy(0, sqrt(N)/sqrt(x_j' x_j))
  for (j in 1:P) {
    real scale_theta = sqrt(N) / sqrt(xjxj[j]);
    lprior += cauchy_lpdf(theta[j] | 0.0, scale_theta);
  }
}

model {
  target += lprior;

  // ----- Build Sigma = sum_j Ktilde_j + sigma^2 I -----
  matrix[N, N] Sigma = rep_matrix(0.0, N, N);

  for (j in 1:P) {
    // Kj from precomputed distances (N x N)
    matrix[N, N] Kj = se_cov_from_D2(D2[j], alpha[j], lambda[j]);
    for (n in 1:N) Kj[n, n] += 1e-12;

    // Efficient: Ktilde = K - PK - KP + PKP, where P = Q R Q'
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
    Ktilde = 0.5 * (Ktilde + Ktilde');  // stabilize symmetry

    Sigma += Ktilde;
  }

  for (n in 1:N) Sigma[n, n] += square(sigma);

  // Mean: mu + X_c * beta1
  vector[N] mu_vec = rep_vector(mu, N) + X_c * beta1;

  {
    matrix[N, N] L = cholesky_decompose(Sigma);
    target += multi_normal_cholesky_lpdf(y_c | mu_vec, L);
  }
}
