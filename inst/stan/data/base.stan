  int<lower=1> T;      // number of time points
  int<lower=1> T_pre;  // number of pre-treatment time points for treated
  int<lower=1> J;      // number of donor series
  vector[T] y;         // treated
  matrix[T, J] Z;      // donors
  real<lower = 0> kappa;          // for weight prior
  real<lower = 0> pr_rate_sigma;  // rate of sigma prior
