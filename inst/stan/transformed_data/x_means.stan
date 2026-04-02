  matrix[N, K] X_y_mean;
  for (i in 1:N) {
    for (k in 1:K) {
      X_y_mean[i, k] = mean(X_y[i, 1:T_pre[i], k]);
    }
  }
  array[N] matrix[K, J] X_z_mean;
  for (k in 1:K) {
    for (j in 1:J) {
      for (i in 1:N) {
        X_z_mean[i, k, j] = mean(X_z[k, 1:T_pre[i], j]);
      }
    }
  }
