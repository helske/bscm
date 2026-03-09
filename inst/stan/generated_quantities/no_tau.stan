  vector[T] synthetic_y = to_vector(normal_rng(synthetic_mean, sigma));
  
  // Bayesian R2
  real R2 = 0;
  {
    real var_res = square(sigma);
    real var_fit = variance(synthetic_mean);
    R2 = var_fit / (var_fit + var_res);
  }
  
  // for leave-one-out cross-validation
  vector[T_pre] log_lik;
  for (t in 1:T_pre) {
    log_lik[t] = normal_lpdf(y[t] | synthetic_mean[t], sigma);
  }
