  // time-varying regression coefficients
  matrix[T, L] gamma = gamma_raw;
  for (l in 1:L) {
    gamma[, l] = cumulative_sum(tau[l] * gamma[, l]);
    gamma[, l] -= mean(gamma[1:T1, l]);
  }
