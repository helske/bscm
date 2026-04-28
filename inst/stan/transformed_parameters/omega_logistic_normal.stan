  array[N] vector[J] omega;
  for (i in 1:N) {
    omega[i] = softmax(kappa * scale * zero_sum_constrain(eta[i]));
  }
