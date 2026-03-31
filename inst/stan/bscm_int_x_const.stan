// Bayesian SCM
// Intercept
// Covariates with constant effects

data {
#include data/base.stan
#include data/alpha.stan
#include data/beta.stan
}
transformed data {
#include transformed_data/base.stan
#include transformed_data/center_z.stan
#include transformed_data/center_x.stan
}
parameters {
#include parameters/base.stan
#include parameters/alpha.stan
#include parameters/beta.stan
}
model {
#include model/base.stan
#include model/alpha.stan
#include model/beta.stan
  for (i in 1:N) {
     matrix[T_pre[i], K] X = X_y_c[i, 1:T_pre[i]];
     for (j in 1:J) {
       X -= X_z_c[j, 1:T_pre[i], ] * omega[i, j];
     }
     y[1:T_pre[i], i] ~ normal_id_glm(
       append_col(X, Z_c[1:T_pre[i]]), a[i], 
       append_row(beta, omega[i]), sigma[i]
       );
  }
}
generated quantities {
  // actual intercept
  vector[N] alpha;
  matrix[T, N] y_mean;
  for (i in 1:N) {
    alpha[i] = a[i] - Z_mean * omega[i] - X_y_mean[i] * beta + 
    dot_product(omega[i], X_z_mean * beta);
    y_mean[, i] = alpha[i] + X_y[i] * beta + Z * omega[i];
    for (j in 1:J) {
      y_mean[, i] -= omega[i, j] * X_z[j] * beta;
    }
  }
#include generated_quantities/base.stan
}
