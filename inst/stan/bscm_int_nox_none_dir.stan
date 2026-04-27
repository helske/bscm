// Bayesian SCM
// Intercept
// No covariates
// Dirichlet prior on omega

data {
#include data/base.stan
#include data/alpha.stan
}
transformed data {
#include transformed_data/omega_dirichlet.stan
#include transformed_data/z_means.stan
}
parameters {
#include parameters/base.stan
#include parameters/omega_dirichlet.stan
#include parameters/a.stan
}
transformed parameters {
#include transformed_parameters/alpha_z.stan
}
model {
#include model/base.stan
#include model/omega_dirichlet.stan
#include model/a.stan
  for (i in 1:N) {
    int Ti = T_pre[i];
    head(y[i], Ti) ~ normal_id_glm(Z[1:Ti], alpha[i], omega[i], sigma[i]);
  }
}
generated quantities {
  matrix[T, N] y_mean;
  for (i in 1:N) {
    y_mean[, i] = alpha[i] + Z * omega[i];
  }
#include generated_quantities/base.stan
}
