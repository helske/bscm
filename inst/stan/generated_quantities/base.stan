  matrix[T, N] y_rep;
  for (i in 1:N) {
    y_rep[, i] = to_vector(normal_rng(y_mean[, i], sigma[i]));
  }
