// Bayesian SCM
// No intercept
// No covariates

data {
#include data/base.stan
}
transformed data {  
#include transformed_data/base.stan
}
parameters {
#include parameters/base.stan
}
transformed parameters {
#include transformed_parameters/base.stan
}
model {
#include model/base.stan
  y[1:T_pre] ~ normal_id_glm(Z[1:T_pre], 0, omega, sigma);
}
generated quantities {
  vector[T] synthetic_mean = Z * omega;
#include generated_quantities/no_tau.stan
#include generated_quantities/base.stan
}
