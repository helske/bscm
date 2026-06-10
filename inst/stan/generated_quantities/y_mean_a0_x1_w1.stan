  {
#include model/X_z_beta.stan
#include model/X_z_gamma.stan
    for (i in 1:N) {
      y_mean[, i] = X_y[i] * beta + rows_dot_product(X_y[i, , tv_idx], gamma) + 
      (Z - X_z_beta - X_z_gamma) * omega[i];
    }
  }
