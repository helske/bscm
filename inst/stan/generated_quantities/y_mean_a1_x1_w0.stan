  {
#include model/cX_z_beta.stan
    for (i in 1:N) {
      y_mean[, i] = a[i] + cX_y[i] * beta + Z_term * omega[i];
    }
  }
