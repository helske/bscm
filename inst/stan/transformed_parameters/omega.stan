  array[N] vector[J] omega;
  if (dirichlet_omega) {
    for (i in 1:N) {
      omega[i] = omega_[i];
    }
  } else {
    for (i in 1:N) {
      omega[i] = softmax(kappa * scale * zero_sum_constrain(eta[i]));
    }
  }
