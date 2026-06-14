  {
    vector[T] mu;
#include model/X_z_beta.stan
#include model/W_z_gamma.stan
    for (i in 1:N) {
      int Ti = T_pre[i];
      mu = X_y[i] * beta + rows_dot_product(W_y[i], gamma) + 
      Z_term * omega[i];
#include generated_quantities/ar1_gq.stan
    }
  }
