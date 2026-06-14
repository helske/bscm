  {
#include model/X_z_beta.stan
#include model/W_z_gamma.stan
    for (i in 1:N) {
      y_mean[, i] = X_y[i] * beta + rows_dot_product(W_y[i], gamma) + 
      Z_term * omega[i];
    }
  }
