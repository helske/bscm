// Bayesian SCM
// Intercept
// No covariates

data {
#include data/base.stan
#include data/alpha.stan
}

transformed data {
#include transformed_data/base.stan
#include transformed_data/center_z.stan
}

parameters {
#include parameters/base.stan
#include parameters/alpha.stan
}
transformed parameters {
#include transformed_parameters/base.stan
}
model {
#include model/base.stan
#include model/alpha.stan
  y[1:T_pre] ~ normal_id_glm(Z_c, a, omega, sigma);
}
generated quantities {
  // actual intercept
  real alpha = a - Z_mean * omega;
  vector[T] synthetic_mean = alpha + Z * omega;
#include generated_quantities/no_tau.stan
#include generated_quantities/base.stan
}
