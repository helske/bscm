  {
    vector[T] mu;
#include model/X_z_beta.stan
    for (i in 1:N) {
      int Ti = T_pre[i];
      mu = alpha[i] + X_y[i] * beta + (Z - X_z_beta) * omega[i];
#include generated_quantities/ar1_gq.stan
    }
  }
