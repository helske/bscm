  int<lower = 1> L;                 // number of time-varying coefficients
  array[L] int<lower = 1> tv_idx;   // columns of X matching tv coefs
  vector<lower = 0>[L] pr_rate_sigma_gamma; // sigma_gamma ~ Gamma(2, .)
  int<lower = 4> D;   // number of spline coefficients + 1
  matrix[T, D - 1] A; // spline projection matrix
