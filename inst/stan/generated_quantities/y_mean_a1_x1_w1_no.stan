  for (i in 1:N) {
    y_mean[, i] = a[i] + cX_y[i] * beta + 
    rows_dot_product(W_y[i], gamma);
  }
