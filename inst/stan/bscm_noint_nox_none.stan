// Bayesian SCM
// No intercept
// No covariates
functions {
#include functions/omega_logistic_normal.stan
}
data {
#include data/base.stan
}
transformed data {
#include transformed_data/omega_logistic_normal.stan
}
parameters {
#include parameters/base.stan
#include parameters/omega_logistic_normal.stan
}
transformed parameters {
#include transformed_parameters/omega_logistic_normal.stan
}
model {
#include model/base.stan
#include model/omega_logistic_normal.stan
  for (i in 1:N) {
    int Ti = T_pre[i];
    head(y[i], Ti) ~ normal_id_glm(Z[1:Ti], 0, omega[i], sigma[i]);
  }
}
generated quantities {
  matrix[T, N] y_mean;
  for (i in 1:N) {
    y_mean[, i] = Z * omega[i];
  }
#include generated_quantities/base.stan
}
