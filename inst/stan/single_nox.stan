// Bayesian SCM with single treated unit and no time-varying covariates
data {
  int<lower=1> T;          // number of time points
  int<lower=1> T_pre;        // Number of pre-treatment time points for treated
  int<lower=1> J;          // number of donor series
  vector[T] y;             // treated
  matrix[T, J] Z;          // donors
  real pr_scale_sigma;
  real pr_mean_intercept; 
  real pr_sd_intercept; 
  real<lower=0> kappa;
}
transformed data {
  // center controls using pretreatment period, levels accounted by intercept
  matrix[T, J] Z_c = Z;
  for (i in 1:J) {
    Z_c[, i] -= mean(Z[1:T_pre, i]);
  }
  // auxiliary stuff
  vector[T - T_pre] tt = linspaced_vector(T - T_pre, 1, T - T_pre);
  vector[J] kappavec = rep_vector(kappa, J);
}
parameters {
  real<lower=0> sigma; // standard deviation of the error term
  real alpha; // intercept
  simplex[J] omega;
}
transformed parameters {
  vector[T] synthetic_mean = alpha + Z_c * omega;
}
model {
  sigma ~ exponential(pr_scale_sigma);
  alpha ~ normal(pr_mean_intercept, pr_sd_intercept);
  omega ~ dirichlet(kappavec);
  y[1:T_pre] ~ normal(synthetic_mean[1:T_pre], sigma);
}
generated quantities {
  vector[T] synthetic_y = to_vector(normal_rng(synthetic_mean, sigma));
  vector[T] effect = y - synthetic_y;
  real avg_effect_pre = mean(head(effect, T_pre));
  real avg_effect_post = mean(tail(effect, T - T_pre));
  vector[T - T_pre] avg_effect_post_cumulative =  cumulative_sum(tail(effect, T - T_pre)) ./ tt;
  vector[T - T_pre] relative_change = tail(synthetic_y, T - T_pre) ./ synthetic_y[T_pre];
  real pre_RMSE = sqrt(mean(head(effect, T_pre)^2));
  real post_RMSE = sqrt(mean(tail(effect, T - T_pre)^2));
  real se = sigma * sqrt(1 + sum(square(omega)));
  real R2 = 0;
  {
    real var_res = square(sigma);
    real var_fit = variance(synthetic_mean);
    R2 = var_fit / (var_fit + var_res);
  }
}
