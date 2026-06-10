if (cv == 0) {
  for (i in 1:N) {
    y_rep[, i] = to_vector(normal_rng(y_mean[, i], sigma[i]));
  }
}
