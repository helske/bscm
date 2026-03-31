// Bayesian SCM
// No intercept
// Covariates with constant effects

data {
#include data/base.stan
#include data/beta.stan
}
transformed data {
#include transformed_data/base.stan
}
parameters {
#include parameters/base.stan
#include parameters/beta.stan
}
model {
#include model/base.stan
#include model/beta.stan
  for (i in 1:N) {
    matrix[T_pre[i], K] X = X_y[i, 1:T_pre[i]];
    for (j in 1:J) {
      X -= X_z[j, 1:T_pre[i], ] * omega[i, j];
    }
    y[1:T_pre[i], i] ~ normal_id_glm(
      append_col(X, Z[1:T_pre[i]]), 0, append_row(beta, omega[i]), sigma[i]
      );
  }
 
}
generated quantities {
  matrix[T, N] y_mean;
  for (i in 1:N) {
    y_mean[, i] = X_y[i] * beta + Z * omega[i];
    for (j in 1:J) {
      y_mean[, i] -= omega[i, j] * X_z[j] * beta;
    }
  }
#include generated_quantities/base.stan
}
