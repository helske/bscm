  {
    vector[T] mu;
    for (i in 1:N) {
      int Ti = T_pre[i];
      mu = alpha[i] + X_y[i] * beta + rows_dot_product(X_y[i, , tv_idx], gamma);
#include generated_quantities/ar1_gq.stan
    }
  }
