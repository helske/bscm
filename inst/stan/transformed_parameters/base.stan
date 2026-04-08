  array[N] vector[J] omega;
  for (i in 1:N) {
    omega[i] = softmax(kappa * eta[i]);
  }
