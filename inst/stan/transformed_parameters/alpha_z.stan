  vector[N] alpha = a;
  for (i in 1:N) {
    alpha[i] -= Z_mean[i] * omega[i];
  }
