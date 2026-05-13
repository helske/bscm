// Regression model, no donors, time-varying coefficients
data {
  int<lower = 1> J;                       // number of donor series
  int<lower = 1> N;                       // number of treated series
  int<lower = 1> T;                       // number of time points
  array[N] int<lower=1> T_pre;            // number of pre-treatment time points
  array[N] vector[T] y;                   // treated
  vector<lower = 0>[N] pr_rate_sigma;     // sigma ~ exponential()
  vector[N] pr_mean_intercept;                 // prior mean of intercept
  vector<lower = 0>[N] pr_sd_intercept;        // prior SD of intercept
  int<lower = 1> K;                // number of covariates
  array[N] matrix[T, K] X_y;       // covariates for treated
  vector[K] pr_mean_beta;          // prior mean of beta
  vector<lower = 0>[K] pr_sd_beta; // prior SD of beta
  int<lower = 1> L;                 // number of time-varying coefficients
  array[L] int<lower = 1> tv_idx;   // columns of X matching tv coefs
  vector<lower = 0>[L] pr_rate_sigma_gamma; // sigma_gamma ~ Gamma(2, .)
}
transformed data {
  matrix[N, K] X_y_mean;
  for (i in 1:N) {
    for (k in 1:K) {
      X_y_mean[i, k] = mean(X_y[i, 1:T_pre[i], k]);
    }
  }
  int T1 = max(T_pre);
  int T2 = T - T1;
}
parameters {
  vector<lower = 0>[N] sigma; // SD of the error term
  vector[N] a; // working intercept
  vector[K] beta; // regression coefficients
  matrix[T, L] gamma_raw;
  vector<lower = 0>[L] sigma_gamma;
}
transformed parameters {
  vector[N] alpha = a - X_y_mean * beta;
  // time-varying regression coefficients
  matrix[T, L] gamma = gamma_raw;
  for (l in 1:L) {
    gamma[, l] = cumulative_sum(sigma_gamma[l] * gamma[, l]);
    gamma[, l] -= mean(gamma[1:T1, l]);
  }
}
model {
  sigma ~ exponential(pr_rate_sigma);
  a ~ normal(pr_mean_intercept, pr_sd_intercept);
  beta ~ normal(pr_mean_beta, pr_sd_beta);
  sigma_gamma ~ gamma(2, pr_rate_sigma_gamma);
  to_vector(gamma_raw) ~ std_normal();
  for (i in 1:N) {
    int Ti = T_pre[i];
    head(y[i], Ti) ~ normal_id_glm(
      X_y[i, 1:Ti], 
      alpha[i] + rows_dot_product(X_y[i, 1:Ti, tv_idx], gamma[1:Ti]), 
      beta, sigma[i]
      );
  }
}
generated quantities {
  matrix[T, N] y_mean;
  for (i in 1:N) {
    y_mean[, i] = alpha[i] + X_y[i] * beta + 
    rows_dot_product(X_y[i, , tv_idx], gamma);
  }
  matrix[T, N] y_rep;
  for (i in 1:N) {
    y_rep[, i] = to_vector(normal_rng(y_mean[, i], sigma[i]));
  }
}
