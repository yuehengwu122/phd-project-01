functions {
  // 1D sine basis on [-L, L], shifted input x in that domain
  // phi_m(x) = sin( m*pi*(x+L)/(2L) ) / sqrt(L),  m = 1..M
  vector phi_sine_1d(vector x, int m, real L) {
    int N = num_elements(x);
    vector[N] out = sin(pi() * m * (x + L) / (2 * L)) / sqrt(L);
    return out;
  }
}

data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;
  vector[N] y;
  int<lower=1> M;
}

transformed data {
  // ---- Center y and X ----
  real y_bar = mean(y);
  vector[N] y_c = y - y_bar;

  row_vector[P] x_bar;
  matrix[N, P] X_c;
  for (p in 1:P) {
    x_bar[p] = mean(col(X, p));
    X_c[, p] = col(X, p) - x_bar[p];
  }

  // ---- dependent prior ingredients: (X'X)^{-1} for centered X ----
  matrix[P, P] XtX = crossprod(X_c);
  matrix[P, P] XtX_ridge = XtX + 1e-10 * diag_matrix(rep_vector(1.0, P));
  matrix[P, P] V_theta_base = inverse_spd(XtX_ridge);
  matrix[P, P] L_V_theta_base = cholesky_decompose(V_theta_base);

  // ---- projector ingredients for each predictor: span{1, x_c} ----
  array[P] matrix[N, 2] Q_tr;
  array[P] matrix[2, 2] R_tr;
  for (p in 1:P) {
    Q_tr[p][, 1] = rep_vector(1.0, N);
    Q_tr[p][, 2] = X_c[, p];
    {
      matrix[2, 2] G = crossprod(Q_tr[p]);   // Q'Q
      R_tr[p] = inverse_spd(G);
    }
  }

  // ---- boundary L per predictor ----
  vector[P] Lp;
  for (p in 1:P) {
    vector[N] xcp = to_vector(X_c[, p]);
    real a = min(xcp);
    real b = max(xcp);
    real half_range = 0.5 * (b - a);
    Lp[p] = fmax(1e-6, 1.5 * half_range);
  }

  // ---- basis PHI, orthogonalized to span{1, x_c} ----
  array[P] matrix[N, M] PHI;
  for (p in 1:P) {
    real L = Lp[p];
    vector[N] xcp = to_vector(X_c[, p]);

    for (m in 1:M)
      PHI[p][, m] = phi_sine_1d(xcp, m, L);

    // PHI <- (I - Q R Q') PHI
    {
      matrix[N, 2] Q = Q_tr[p];
      matrix[2, 2] R = R_tr[p];
      PHI[p] = PHI[p] - Q * (R * (Q' * PHI[p]));
    }
  }

  // ---- frequencies omega_m(p) = m*pi/(2*Lp[p]) ----
  array[P] vector[M] omega_pm;
  for (p in 1:P) {
    real L = Lp[p];
    for (m in 1:M)
      omega_pm[p][m] = m * pi() / (2 * L);
  }
}

parameters {
  // dependent (mixture-of-g / JZS-style) prior scale for theta
  // real<lower=0> g_theta;
  real log_g_theta;
  vector[P] z_theta;         // theta = sqrt(g_theta) * L * z_theta

  real mu;                   // intercept on centered y
  real<lower=0> sigma;       // noise SD

  vector<lower=0>[P] delta;  // standardized GP amplitude alpha/sigma
  vector<lower=0>[P] lambda; // length-scale

  array[P] vector[M] z_basis; // standard normal weights for HSGP
}

transformed parameters {
  // theta | g_theta ~ N(0, g_theta * (X_c'X_c)^{-1})
  // vector[P] theta = sqrt(g_theta) * (L_V_theta_base * z_theta);

  real<lower=0> g_theta = exp(log_g_theta);
  vector[P] theta = exp(0.5 * log_g_theta) * (L_V_theta_base * z_theta);

  // slopes on y-scale
  vector[P] beta1 = sigma * theta;

  // prior buffer (useful for bridge sampling)
  real lprior = 0;

  // match your previous priors
  for (p in 1:P)
    lprior += student_t_lpdf(delta[p] | 4, 0, 2.7) - student_t_lccdf(0 | 4, 0, 2.7);

  for (p in 1:P)
    lprior += inv_gamma_lpdf(lambda[p] | 3, 0.8);

  lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(mu | 3, 0, 10);

  // same g prior as your exact model (2)
  // lprior += inv_gamma_lpdf(g_theta | 0.5, 0.5 * N);
  lprior += normal_lpdf(log(g_theta) | log(N), 0.5);
}

model {
  target += lprior;

  z_theta ~ std_normal();
  for (p in 1:P) z_basis[p] ~ std_normal();

  // mean = linear + HSGP nonlinear
  vector[N] mu_vec = rep_vector(mu, N) + X_c * beta1;

  for (p in 1:P) {
    // HSGP spectral weights for Matérn-3/2 (same structure as your (1))
    real a = delta[p] * sigma;   // alpha = delta*sigma
    real ell = lambda[p];

    vector[M] w;
    for (m in 1:M) {
      real om = omega_pm[p][m];
      real denom = square(1 + square(ell * om));  // (1 + (ell*ω)^2)^2
      real s = square(a) * ell / denom;           // scaling constant absorbed by delta
      w[m] = sqrt(s) * z_basis[p][m];
    }

    mu_vec += PHI[p] * w;
  }

  y_c ~ normal(mu_vec, sigma);
}
