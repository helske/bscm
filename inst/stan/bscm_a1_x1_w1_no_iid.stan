// Regression model, no donors, time-varying coefficients
data {
#include data/base.stan
#include data/alpha.stan
#include data/beta.stan
#include data/gamma.stan
}
transformed data {
  matrix[N, K] X_y_mean;
  for (i in 1:N) {
    for (k in 1:K) {
      X_y_mean[i, k] = mean(X_y[i, 1:T_pre[i], k]);
    }
  }
  int T1 = max(T_pre);
  int T2 = T - T1;
}
parameters {
#include parameters/base.stan
#include parameters/a.stan
#include parameters/beta.stan
#include parameters/gamma.stan
}
transformed parameters {
#include transformed_parameters/gamma.stan
  vector[N] alpha = a - X_y_mean * beta;
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
#include generated_quantities/y_mean_a1_x1_w1_no.stan
#include generated_quantities/y_rep_iid.stan
}
