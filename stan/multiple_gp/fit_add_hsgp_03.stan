functions {
  // 1D sine basis on [-L, L], shifted input x in that domain
  // phi_m(x) = sin( m*pi*(x+L)/(2L) ) / sqrt(L)
  // m = 1..M
  vector phi_sine_1d(vector x, int m, real L) {
    int N = num_elements(x);
    vector[N] out;
    out = sin(pi() * m * (x + L) / (2 * L)) / sqrt(L);
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
  // center y and X
  real y_bar = mean(y);
  vector[N] y_c = y - y_bar;

  row_vector[P] x_bar;
  matrix[N, P] X_c;
  for (p in 1:P) {
    x_bar[p] = mean(col(X, p));
    X_c[, p] = col(X, p) - x_bar[p];
  }

  // projector ingredients for each predictor: span{1, x_c}
  array[P] matrix[N, 2] Q_tr;
  array[P] matrix[2, 2] R_tr;
  vector[P] xjxj;

  for (p in 1:P) {
    Q_tr[p][, 1] = rep_vector(1.0, N);
    Q_tr[p][, 2] = X_c[, p];
    {
      matrix[2, 2] G = crossprod(Q_tr[p]);
      R_tr[p] = inverse_spd(G);
    }
    xjxj[p] = dot_self(X_c[, p]);
  }

  // boundary L per predictor
  vector[P] Lp;
  for (p in 1:P) {
    vector[N] xcp = to_vector(X_c[, p]);
    real a = min(xcp);
    real b = max(xcp);
    real half_range = 0.5 * (b - a);
    Lp[p] = fmax(1e-6, 1.5 * half_range);
  }

  // basis PHI, orthogonalized to span{1, x_c}
  array[P] matrix[N, M] PHI;
  for (p in 1:P) {
    real L = Lp[p];
    vector[N] xcp = to_vector(X_c[, p]);

    for (m in 1:M) {
      PHI[p][, m] = phi_sine_1d(xcp, m, L);
    }

    // PHI <- (I - Q R Q') PHI
    {
      matrix[N, 2] Q = Q_tr[p];
      matrix[2, 2] R = R_tr[p];
      PHI[p] = PHI[p] - Q * (R * (Q' * PHI[p]));
    }
  }

  // Precompute frequencies omega_m(p) = m*pi/(2*Lp[p]) for all p,m
  array[P] vector[M] omega_pm;
  for (p in 1:P) {
    real L = Lp[p];
    for (m in 1:M)
      omega_pm[p][m] = m * pi() / (2 * L);
  }
}

parameters {
  real mu;                 // intercept on centered y
  vector[P] theta;         // standardized slopes beta/sigma
  real<lower=0> sigma;     // noise SD

  vector<lower=0>[P] delta;   // standardized GP amplitude alpha/sigma
  vector<lower=0>[P] lambda;  // length-scale (Matérn parameterization)

  array[P] vector[M] z_basis; // standard normal weights
}

transformed parameters {
  vector[P] beta1 = sigma * theta;        // slopes on y-scale
  real lprior = 0;

  // priors (match your previous)
  for (p in 1:P)
    lprior += student_t_lpdf(delta[p] | 4, 0, 2.7) - student_t_lccdf(0 | 4, 0, 2.7);

  for (p in 1:P)
    lprior += inv_gamma_lpdf(lambda[p] | 3, 1.5);

  lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);

  lprior += student_t_lpdf(mu | 3, 0, 10);

  for (p in 1:P) {
    real scale_theta = sqrt(N) / sqrt(xjxj[p]);
    lprior += cauchy_lpdf(theta[p] | 0.0, scale_theta);
  }
}

model {
  target += lprior;

  for (p in 1:P) z_basis[p] ~ std_normal();

  // mean = linear + HSGP nonlinear
  vector[N] mu_vec = rep_vector(mu, N) + X_c * beta1;

  for (p in 1:P) {
    real a = delta[p] * sigma;   // alpha = delta*sigma
    real ell = lambda[p];        // SE length-scale

    // --- SE spectrum shape (no amplitude), then normalize so sum(s_norm)=1 ---
    vector[M] s_raw;
    for (m in 1:M) {
      real om = omega_pm[p][m];
      s_raw[m] = ell * exp(-0.5 * square(ell * om));  // shape-only spectrum
    }
    real s_sum = sum(s_raw) + 1e-12;  // numerical guard

    vector[M] w;
    for (m in 1:M) {
      real s = square(a) * (s_raw[m] / s_sum);  // amplitude only via a^2
      w[m] = sqrt(s) * z_basis[p][m];
    }

    mu_vec += PHI[p] * w;
  }

  y_c ~ normal(mu_vec, sigma);
}
