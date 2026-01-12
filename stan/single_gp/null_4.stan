data {
  int<lower=1> N;
  array[N] real x;        // 1D inputs
  vector[N] y;            // response
}

transformed data {
  // Build X = [1, x] for the linear (ridge) term without changing the data block
  matrix[N, 2] X;
  for (n in 1:N) {
    X[n, 1] = 1;
    X[n, 2] = x[n];
  }
}

parameters {
  real<lower=0> tau2;
  real<lower=0> sigma;    // noise SD
}

transformed parameters {

  // Buffer fully normalized priors for bridge sampling
  real lprior = 0;
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)  - student_t_lccdf(0 | 3, 0, 2.5);
  // tau2 ~ IG(0.5, 0.5*s^2), s is scale of Cauchy distribution, fixed at 0.7 here
  lprior += inv_gamma_lpdf(tau2 | 0.5, 0.5 * 0.7 * 0.7); 
}

model {
  target += lprior;

  // Σ_y = σ^2 τ^2 X (X'X)^{-1} X' + K(x; α, λ) + σ^2 I

  matrix[2, 2] XX_inv = inverse_spd(crossprod(X)); // (X'X)^{-1}
  matrix[N, N] Sigma = square(sigma) * tau2 * quad_form_sym(XX_inv, X');
  for (n in 1:N) Sigma[n, n] += square(sigma);

  matrix[N, N] L = cholesky_decompose(Sigma);
  target += multi_normal_cholesky_lpdf(y | rep_vector(0.0, N), L);
}
