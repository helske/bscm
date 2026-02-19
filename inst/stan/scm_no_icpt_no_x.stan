// Bayesian SCM with single treated unit and no intercept or covariates
data {
  int<lower=1> T;          // number of time points
  int<lower=1> T_pre;        // Number of pre-treatment time points for treated
  int<lower=1> J;          // number of donor series
  vector[T] y;             // treated
  matrix[T, J] Z;          // donors
  real pr_rate_sigma;
  real<lower = 0> kappa;
}
transformed data {  
  // auxiliary stuff
  vector[T - T_pre] tt = linspaced_vector(T - T_pre, 1, T - T_pre);
}
parameters {
  real<lower = 0> sigma; // standard deviation of the error term
  vector<lower = 0>[J] omega_raw; // Gamma(kappa, 1);
}
transformed parameters {
  // to simplex, omega ~ Dirichlet(kappa);
  vector[J] omega = omega_raw / sum(omega_raw);
  vector[T] synthetic_mean = Z * omega;
}
model {
  sigma ~ exponential(pr_rate_sigma);
  omega_raw ~ gamma(kappa, 1);
  y[1:T_pre] ~ normal(synthetic_mean[1:T_pre], sigma);
}
generated quantities {
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
