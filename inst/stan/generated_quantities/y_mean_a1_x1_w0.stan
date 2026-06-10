  {
#include model/X_z_beta.stan
    for (i in 1:N) {
      y_mean[, i] = alpha[i] + X_y[i] * beta + (Z - X_z_beta) * omega[i];
    }
  }
