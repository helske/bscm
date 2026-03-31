  int<lower=1> J;      // number of donor series
  int<lower=1> N;      // number of treated series
  int<lower=1> T;      // number of time points
  array[N] int<lower=1> T_pre;  // number of pre-treatment time points
  matrix[T, N] y;      // treated
  matrix[T, J] Z;      // donors
  real<lower = 0> kappa;               // for weight prior
  vector<lower = 0>[N] pr_rate_sigma;  // rate of sigma prior
