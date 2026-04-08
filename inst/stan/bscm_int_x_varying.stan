// Bayesian SCM
// Intercept
// Covariates with varying effects

data {
#include data/base.stan
#include data/alpha.stan
#include data/beta.stan
#include data/gamma.stan
}
transformed data {
#include transformed_data/base.stan
#include transformed_data/z_means.stan
#include transformed_data/x_means.stan
#include transformed_data/gamma.stan
}
parameters {
#include parameters/base.stan
#include parameters/a.stan
#include parameters/beta.stan
#include parameters/gamma.stan
}
transformed parameters {
#include transformed_parameters/base.stan
#include transformed_parameters/alpha_z.stan
#include transformed_parameters/alpha_x.stan
#include transformed_parameters/gamma.stan
}
model {
#include model/base.stan
#include model/a.stan
#include model/beta.stan
#include model/gamma.stan
  {
    matrix[T, K] X;
    for (i in 1:N) {
      int Ti = T_pre[i];
      X = X_y[i];
      for (k in 1:K) {
        X[, k] -= X_z[k] * omega[i];
      }
      head(y[i], Ti) ~ normal_id_glm(
        Z[1:Ti], alpha[i] + X[1:Ti] * beta +
        rows_dot_product(X[1:Ti, tv_idx], gamma[1:Ti]),
        omega[i], sigma[i]
        );
    }
  }
}
generated quantities {
  matrix[T, N] y_mean;
  for (i in 1:N) {
    y_mean[, i] = alpha[i] + X_y[i] * beta + Z * omega[i] +
    rows_dot_product(X_y[i, , tv_idx], gamma);
    for (k in 1:K) {
      y_mean[, i] -= X_z[k] * omega[i] * beta[k];
    }
    for (l in 1:L) {
      y_mean[, i] -= X_z[tv_idx[l]] * omega[i] .* gamma[, l];
    }
  }
#include generated_quantities/base.stan
}
