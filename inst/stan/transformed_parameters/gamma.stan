  // time-varying regression coefficients
  matrix[T, L] gamma;
  for (l in 1:L) {
    gamma[, l] = A * (sigma_gamma[l] * xi[, l]);
  }
