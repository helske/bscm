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
  matrix[N, K] X_y_mean;
  for (i in 1:N) {
    for (k in 1:K) {
      X_y_mean[i, k] = mean(X_y[i, 1:T_pre[i], k]);
    }
  }
}
parameters {
#include parameters/base.stan
#include parameters/a.stan
#include parameters/beta.stan
#include parameters/rho.stan
}
transformed parameters {
  vector[N] alpha = a - X_y_mean * beta;
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
#include generated_quantities/gq_ar1_a1_x1_w0_no.stan
}
