if (likelihood) {
  vector[T] mu;
#include model/X_z_beta.stan
#include model/X_z_gamma.stan
  for (i in 1:N) {
    int Ti = T_pre[i];
    mu = X_y[i] * beta + rows_dot_product(X_y[i, , tv_idx], gamma) + 
    (Z - X_z_beta - X_z_gamma) * omega[i];
    head(y[i], Ti) ~ normal(head(mu, Ti), sigma[i]);
  }
}
