  alpha -= X_y_mean * beta;
  {
    vector[J] xzb = X_z_mean * beta;
    for (i in 1:N) {
      alpha[i] += dot_product(xzb, omega[i]);
    }
  }
