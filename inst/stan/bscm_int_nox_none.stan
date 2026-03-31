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
model {
#include model/base.stan
#include model/alpha.stan
  for (i in 1:N) {
    y[1:T_pre[i], i] ~ normal_id_glm(Z_c[1:T_pre[i]], a[i], omega[i],  sigma[i]);
  }
}
generated quantities {
  // actual intercept
  vector[N] alpha;
  matrix[T, N] y_mean;
  for (i in 1:N) {
    alpha[i] = a[i] - Z_mean * omega[i];
    y_mean[, i] = alpha[i] + Z * omega[i];
  }
#include generated_quantities/base.stan
}
