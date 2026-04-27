// Bayesian SCM
// No intercept
// Covariates with constant effects
// Dirichlet prior on omega

data {
#include data/base.stan
#include data/beta.stan
}
transformed data {
#include transformed_data/omega_dirichlet.stan
}
parameters {
#include parameters/base.stan
#include parameters/omega_dirichlet.stan
#include parameters/beta.stan
}
model {
#include model/base.stan
#include model/omega_dirichlet.stan
#include model/beta.stan
  {
    matrix[T, K] X;
    for (i in 1:N) {
      int Ti = T_pre[i];
      X = X_y[i];
      for (k in 1:K) {
        X[, k] -= X_z[k] * omega[i];
      }
      head(y[i], Ti) ~ normal_id_glm(
        Z[1:Ti], X[1:Ti] * beta, omega[i], sigma[i]
        );
    }
  }
}
generated quantities {
  matrix[T, N] y_mean;
  for (i in 1:N) {
    y_mean[, i] = X_y[i] * beta + Z * omega[i];
    for (k in 1:K) {
      y_mean[, i] -= X_z[k] * omega[i] * beta[k];
    }
  }
#include generated_quantities/base.stan
}
