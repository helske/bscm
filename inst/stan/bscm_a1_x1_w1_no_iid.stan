// Regression model, no donors, time-varying coefficients
data {
#include data/base.stan
#include data/alpha.stan
#include data/beta.stan
#include data/gamma.stan
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
  array[N] matrix[T, L] W_y;
  for (i in 1:N) {
    W_y[i] = X_y[i, , tv_idx];
  }
}
parameters {
#include parameters/base.stan
#include parameters/a.stan
#include parameters/beta.stan
#include parameters/gamma.stan
}
transformed parameters {
#include transformed_parameters/gamma.stan
}
model {
#include model/base.stan
#include model/a.stan
#include model/beta.stan
#include model/gamma.stan
#include model/likelihood_a1_x1_w1_no_iid.stan
}
generated quantities {
#include generated_quantities/base.stan
  vector[N] alpha = a - X_y_mean * beta;
#include generated_quantities/y_mean_a1_x1_w1_no.stan
#include generated_quantities/y_rep_iid.stan
}
