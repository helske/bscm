  {
    vector[T] mu;
#include model/cX_z_beta.stan
    for (i in 1:N) {
      int Ti = T_pre[i];
      mu = a[i] + cX_y[i] * beta + Z_term * omega[i];
#include generated_quantities/ar1_gq.stan
    }
  }
