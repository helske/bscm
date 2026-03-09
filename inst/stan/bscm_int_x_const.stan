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
transformed parameters {
#include transformed_parameters/base.stan
}
model {
#include model/base.stan
#include model/alpha.stan
#include model/beta.stan
  {
     matrix[T_pre, K] X = X_y_c;
     for (j in 1:J) {
       X -= X_z_c[j] * omega[j];
     }
     y[1:T_pre] ~ normal_id_glm(
       append_col(X, Z_c), a, append_row(beta, omega), sigma
       );
  }
}
generated quantities {
  real alpha = a - Z_mean * omega - X_y_mean * beta + 
  dot_product(omega, X_z_mean * beta);
  vector[T] synthetic_mean = alpha + X_y * beta + Z * omega;
  for (j in 1:J) {
    synthetic_mean -= omega[j] * X_z[j] * beta;
  }
#include generated_quantities/no_tau.stan
#include generated_quantities/base.stan
}
