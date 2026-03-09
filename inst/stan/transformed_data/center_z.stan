  // center donors using the pre-treatment period
  // uncertainty of alpha_z is absorbed to a_y
  matrix[T_pre, J] Z_c = Z[1:T_pre];
  row_vector[J] Z_mean;
  for (j in 1:J) {
    Z_mean[j] = mean(Z_c[, j]);
    Z_c[, j] -= Z_mean[j];
  }
