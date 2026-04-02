  alpha -= X_y_mean * beta;
  for (i in 1:N) {
    alpha[i] += dot_product(X_z_mean[i] * omega[i], beta);
  }
