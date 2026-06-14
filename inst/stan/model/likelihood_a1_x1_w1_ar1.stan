  if (likelihood) {
    vector[T] mu;
#include model/cX_z_beta.stan
#include model/W_z_gamma.stan
    for (i in 1:N) {
      int Ti = T_pre[i];
      mu = a[i] + cX_y[i] * beta + rows_dot_product(W_y[i], gamma) + 
      Z_term * omega[i];
      mu[2:Ti] += rho[i] * head(y[i] - mu, Ti - 1);
      head(y[i], Ti) ~ normal(head(mu, Ti), sigma[i]);
    }
  }
