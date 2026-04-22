  matrix[T, N] y_rep;
  matrix[T, N] effect;
  for (i in 1:N) {
    y_rep[, i] = to_vector(normal_rng(y_mean[, i], sigma[i]));
    effect[, i] = y[i] - y_rep[, i];
  }
