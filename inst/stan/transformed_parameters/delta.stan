  vector[T] delta;
  {
    matrix[T, K] X = X_y;
    for (j in 1:J) {
      X -= X_z[j] * omega[j];
    }
    delta = rows_dot_self(diag_post_multiply(X, inv_sd_x));
  }
