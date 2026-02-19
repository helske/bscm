// Bayesian SCM with single treated unit and time-varying intercept
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
  vector<lower = 0>[J] pr_rate_sigma_z;
  vector[J] pr_mean_intercept_z;
  vector<lower = 0>[J] pr_sd_intercept_z;
  real<lower = 0> pr_sd_sigma_delta;
  int<lower = 0> D; // degrees of freedom of the spline
  matrix[T, D] spline_matrix; // transformed spline basis matrix
}
transformed data {
  // auxiliary stuff
  vector[T - T_pre] tt = linspaced_vector(T - T_pre, 1, T - T_pre);
  
}
parameters {
  real<lower = 0> sigma; // standard deviation of the error term
  vector<lower = 0>[J] sigma_z; // SD of the error terms of donor equations
  vector<lower = 0>[J] omega_raw; // Gamma(kappa, 1);
  real alpha; // intercept
  vector[J] alpha_z; // intercepts for donors
  vector[D] delta; // spline coefficients
  real<lower = 0> sigma_delta;
}
transformed parameters {
  // to simplex, omega ~ Dirichlet(kappa);
  vector[J] omega = omega_raw / sum(omega_raw);
  //vector[T] mu = spline_matrix * cumulative_sum(sigma_delta * delta);
  vector[T] mu = spline_matrix * sigma_delta * delta;
}
model {
  sigma_delta ~ normal(0, pr_sd_sigma_delta);
  //delta ~ std_normal();
  delta[1] ~ normal(0, 1);
  delta[2:D] ~ normal(delta[1:(D - 1)], 1);
  alpha_z ~ normal(pr_mean_intercept_z, pr_sd_intercept_z);
  sigma_z ~ exponential(pr_rate_sigma_z);
  for (j in 1:J) {
    Z[, j] ~ normal(mu + alpha_z[j], sigma_z[j]);
  }
  omega_raw ~ gamma(kappa, 1);
  alpha ~ normal(pr_mean_intercept, pr_sd_intercept);
  sigma ~ exponential(pr_rate_sigma);
  {
    // mu term cancels out
    vector[T_pre] Z_omega = Z[1:T_pre, ] * omega;
    for (j in 1:J) {
      Z_omega -= alpha_z[j] * omega[j];
    }
    y[1:T_pre] ~ normal(alpha + Z_omega, sigma);
  }
}
generated quantities {
  vector[T] synthetic_mean = alpha + Z * omega;
  for (j in 1:J) {
    synthetic_mean -= alpha_z[j] * omega[j];
  }
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
