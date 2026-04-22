  int<lower = 1> J;                       // number of donor series
  int<lower = 1> N;                       // number of treated series
  int<lower = 1> T;                       // number of time points
  array[N] int<lower=1> T_pre;            // number of pre-treatment time points
  array[N] vector[T] y;                   // treated
  matrix[T, J] Z;                         // donors
  vector<lower = 0>[N] pr_rate_sigma;     // sigma ~ exponential()
  real<lower = 0> kappa;
