  // center donors using the pre-treatment period
  // uncertainty of alpha_z is absorbed to a_y
  matrix[T, J] Z_c;
  row_vector[J] Z_mean;
  {
    int min_T_pre = min(T_pre);
    for (j in 1:J) {
      Z_mean[j] = mean(Z[1:min_T_pre, j]);
      Z_c[, j] = Z[, j] - Z_mean[j];
    }
  }
