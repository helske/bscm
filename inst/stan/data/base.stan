  int<lower = 1> J;                       // number of donor series
  int<lower = 1> N;                       // number of treated series
  int<lower = 1> T;                       // number of time points
  array[N] int<lower=1> T_pre;            // number of pre-treatment time points
  array[N] vector[T] y;                   // treated
  matrix[T, J] Z;                         // donors
  real<lower = 0> kappa;                  // hyperparameter of prior on omega
  vector<lower = 0>[N] pr_rate_sigma;     // sigma ~ exponential()
  int<lower = 0, upper = 1> cv;           // skip computing y_rep if 1
  int<lower = 0, upper = 1> likelihood;   // if 0, sample from prior predictive
  int<lower = 0, upper = 1> dirichlet_omega;  // dirichlet or logistic normal
