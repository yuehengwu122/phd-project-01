data {
  int<lower=1> N;
  array[N] real x;
  real<lower=0> lambda;
  real<lower=0> sigma;
  real<lower=0> delta;
}

transformed data {
  real alpha = delta * sigma;

  matrix[N, N] cov = gp_exp_quad_cov(x, alpha, lambda)
                     + diag_matrix(rep_vector(1e-10, N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
}

generated quantities {
  vector[N] f = multi_normal_cholesky_rng(rep_vector(0, N), L_cov);
  array[N] real y = normal_rng(f, sigma);
}
