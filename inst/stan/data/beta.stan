  int<lower = 1> K;                // number of covariates
  array[N] matrix[T, K] X_y;       // covariates for treated
  array[J] matrix[T, K] X_z;       // covariates for donors
  vector[K] pr_mean_beta;          // prior mean of beta
  vector<lower = 0>[K] pr_sd_beta; // prior SD of beta
