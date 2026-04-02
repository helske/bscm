  int<lower = 1> L;                 // number of time-varying coefficients
  array[L] int<lower = 1> tv_idx;   // columns of X matching tv coefs
  vector<lower = 0>[L] pr_rate_tau; // rate for the prior of tau
