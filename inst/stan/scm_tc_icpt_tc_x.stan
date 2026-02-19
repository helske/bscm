// Bayesian SCM with single treated unit and intercept + covariates
data {
  int<lower = 1> T;          // number of time points
  int<lower = 1> T_pre;      // Number of pre-treatment time points for treated
  int<lower = 1> J;          // number of donor series
  vector[T] y;               // treated
  matrix[T, J] Z;            // donors
  real<lower = 0> pr_rate_sigma;
  real pr_mean_intercept; 
  real<lower = 0> pr_sd_intercept; 
  real<lower = 0> kappa;
  int<lower = 1> K;          // number of covariates
  array[J] matrix[T, K] X_z; // covariates for donors
  matrix[T, K] X_y;          // covariates for treated
  vector[K] pr_mean_coef; 
  vector<lower = 0>[K] pr_sd_coef; 
  vector<lower = 0>[J] pr_rate_sigma_z;
  vector[J] pr_mean_intercept_z;
  vector<lower = 0>[J] pr_sd_intercept_z;
}
transformed data {
  // center covariates 
  matrix[T_pre, K] X_y_c = X_y[1:T_pre,];
  row_vector[K] X_y_mean;
  for (k in 1:K) {
    X_y_mean[k] = mean(X_y_c[, k]);
    X_y_c[, k] -= X_y_mean[k];
  }
  array[J] matrix[T, K] X_z_c;
  array[J] row_vector[K] X_z_mean;
  for (j in 1:J) {
    for (k in 1:K) {
      X_z_mean[j, k] = mean(X_z[j, , k]);
      X_z_c[j, , k] = X_z[j, , k] - X_z_mean[j, k];
    }
  }
  // auxiliary stuff
  vector[T - T_pre] tt = linspaced_vector(T - T_pre, 1, T - T_pre);
}
parameters {
  real<lower = 0> sigma; // standard deviation of the error term
  real a; // intercept
  vector<lower = 0>[J] omega_raw; // Gamma(kappa, 1);
  vector[K] beta; // regression coefficients
  vector[J] a_z;
  vector<lower = 0>[J] sigma_z; // SD of the error terms of donor equations
}
transformed parameters {
  // to simplex, omega ~ Dirichlet(kappa);
  vector[J] omega = omega_raw / sum(omega_raw);
}
model {
  beta ~ normal(pr_mean_coef, pr_sd_coef);
  a_z ~ normal(pr_mean_intercept_z, pr_sd_intercept_z);
  sigma_z ~ exponential(pr_rate_sigma_z);
  for (j in 1:J) {
    Z[, j] ~ normal_id_glm(X_z_c[j], a_z[j], beta, sigma_z[j]);
  }
  omega_raw ~ gamma(kappa, 1);
  a ~ normal(pr_mean_intercept, pr_sd_intercept);
  sigma ~ exponential(pr_rate_sigma);
  {
    vector[T_pre] Z_omega = Z[1:T_pre, ] * omega;
    for (j in 1:J) {
      Z_omega -= (a_z[j] + X_z_c[j, 1:T_pre, ] * beta) * omega[j];
    }
    y[1:T_pre] ~ normal_id_glm(X_y_c, a + Z_omega, beta, sigma);
  }
}
generated quantities {
  real alpha = a - X_y_mean * beta;
  vector[T] synthetic_mean = alpha + X_y * beta + Z * omega;
  for (j in 1:J) {
    synthetic_mean -= (a_z[j] + X_z_c[j] * beta) * omega[j];
  }
  vector[J] alpha_z;
  for (j in 1:J) alpha_z[j] = a_z[j] - X_z_mean[j] * beta;
  vector[T] synthetic_y = to_vector(normal_rng(synthetic_mean, sigma));
  vector[T] effect = y - synthetic_y;
  real avg_effect_pre = mean(head(effect, T_pre));
  real avg_effect_post = mean(tail(effect, T - T_pre));
  vector[T - T_pre] avg_effect_post_cumulative = cumulative_sum(tail(effect, T - T_pre)) ./ tt;
  vector[T - T_pre] relative_change = tail(synthetic_y, T - T_pre) ./ synthetic_y[T_pre];
  real RMSE_pre = sqrt(mean(head(effect, T_pre)^2));
  real RMSE_post = sqrt(mean(tail(effect, T - T_pre)^2));
  real RMSE_ratio = RMSE_post / RMSE_pre;
  real effective_donors = inv(sum(square(omega)));
  real R2 = 0;
  {
    real var_res = square(sigma);
    real var_fit = variance(synthetic_mean);
    R2 = var_fit / (var_fit + var_res);
  }
  vector[T_pre] log_lik;
  for (t in 1:T_pre) {
    log_lik[t] = normal_lpdf(y[t] | synthetic_mean[t], sigma);
  }
}
