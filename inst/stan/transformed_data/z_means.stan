  // centered donors
  row_vector[J] Z_mean;
  matrix[T, J] cZ;
  for (j in 1:J) {
    Z_mean[j] = mean(Z[1:T0, j]);
    cZ[, j] = Z[, j] - Z_mean[j];
  }
