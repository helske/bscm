  matrix[T, J] X_z_gamma;
  for (j in 1:J) {
    X_z_gamma[, j] = rows_dot_product(X_z[j][, tv_idx], gamma);
  }
