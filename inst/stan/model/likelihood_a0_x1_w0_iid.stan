  if (likelihood) {
    vector[T] mu;
#include model/X_z_beta.stan
    for (i in 1:N) {
      int Ti = T_pre[i];
      mu = X_y[i] * beta + Z_term * omega[i];
      head(y[i], Ti) ~ normal(head(mu, Ti), sigma[i]);
    }
  }
