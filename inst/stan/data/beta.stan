  int<lower = 1> K;                // number of covariates
  array[N] matrix[T, K] X_y;       // covariates for treated
  array[J] matrix[T, K] X_z;       // covariates for donors
  vector[K] pr_mean_coef;          // prior mean of beta
  vector<lower = 0>[K] pr_sd_coef; // prior SD of beta
  vector<lower = 0>[K] inv_sd_x;   // 1/SD(X_k) over units in t = 1, ..., T_pre
