  matrix[N, K] X_y_mean;
  matrix[J, K] X_z_mean;
  
  array[N] matrix[T, K] cX_y;
  array[J] matrix[T, K] cX_z;
  
  for (i in 1:N) {
    for (k in 1:K) {
      X_y_mean[i, k] = mean(X_y[i, 1:T0, k]);
      cX_y[i, , k] = X_y[i, , k] - X_y_mean[i, k];
    }
  }
  for (j in 1:J) {
    for (k in 1:K) {
      X_z_mean[j, k] = mean(X_z[j, 1:T0, k]);
      cX_z[j, , k] = X_z[j, , k] - X_z_mean[j, k];
    }
  }
