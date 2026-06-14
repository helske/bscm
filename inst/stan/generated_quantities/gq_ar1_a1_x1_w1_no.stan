  {
    vector[T] mu;
    for (i in 1:N) {
      int Ti = T_pre[i];
      mu = a[i] + cX_y[i] * beta + rows_dot_product(W_y[i], gamma);
#include generated_quantities/ar1_gq.stan
    }
  }
