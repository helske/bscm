  for (i in 1:N) {
    y_mean[, i] = alpha[i] + X_y[i] * beta + 
    rows_dot_product(X_y[i, , tv_idx], gamma);
  }
