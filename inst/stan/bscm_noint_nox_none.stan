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
model {
#include model/base.stan
  for (i in 1:N) {
    y[1:T_pre[i], i] ~ normal_id_glm(Z[1:T_pre[i]], 0, omega[i], sigma[i]);
  }
}
generated quantities {
  matrix[T, N] y_mean;
  for (i in 1:N) {
    y_mean[, i] = Z * omega[i];
  }
#include generated_quantities/base.stan
}
