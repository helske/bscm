// Bayesian SCM
// No intercept
// No covariates
// Dirichlet prior on omega

data {
#include data/base.stan
}
transformed data {
#include transformed_data/omega_dirichlet.stan
}
parameters {
#include parameters/base.stan
#include parameters/omega_dirichlet.stan
}
model {
#include model/base.stan
#include model/omega_dirichlet.stan
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
