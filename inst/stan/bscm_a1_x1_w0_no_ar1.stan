// Intercept
// Covariates with constant effects
// No donors
// AR(1) residuals

data {
#include data/base.stan
#include data/alpha.stan
#include data/beta.stan
}
transformed data {
#include transformed_data/base.stan
  matrix[N, K] X_y_mean;
  array[N] matrix[T, K] cX_y;
  
  for (i in 1:N) {
    for (k in 1:K) {
      X_y_mean[i, k] = mean(X_y[i, 1:T0, k]);
      cX_y[i, , k] = X_y[i, , k] - X_y_mean[i, k];
    }
  }
}
parameters {
#include parameters/base.stan
#include parameters/a.stan
#include parameters/beta.stan
#include parameters/rho.stan
}
model {
#include model/base.stan
#include model/a.stan
#include model/beta.stan
#include model/rho.stan
#include model/likelihood_a1_x1_w0_no_ar1.stan
}
generated quantities {
#include generated_quantities/base.stan
  vector[N] alpha = a - X_y_mean * beta;
#include generated_quantities/gq_ar1_a1_x1_w0_no.stan
}
