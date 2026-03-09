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
transformed parameters {
#include transformed_parameters/base.stan
}
model {
#include model/base.stan
#include model/beta.stan
  {
    matrix[T_pre, K] X = X_y[1:T_pre];
    for (j in 1:J) {
      X -= X_z[j, 1:T_pre, ] * omega[j];
    }
    y[1:T_pre] ~ normal_id_glm(
      append_col(X, Z[1:T_pre]), 0, append_row(beta, omega), sigma
      );
  }
 
}
generated quantities {
  vector[T] synthetic_mean = X_y * beta + Z * omega;
  for (j in 1:J) {
    synthetic_mean -= omega[j] * X_z[j] * beta;
  }
#include generated_quantities/no_tau.stan
#include generated_quantities/base.stan
#include transformed_parameters/delta.stan

}
