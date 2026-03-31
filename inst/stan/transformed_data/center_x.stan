  // center covariates using pre-treatment period
  // working intercept a_y absobs the shifts
  
  // centered covariates for treated
  array[N] matrix[T, K] X_y_c;
  matrix[N, K] X_y_mean;
  for (i in 1:N) {
    for (k in 1:K) {
      X_y_mean[i, k] = mean(X_y[i, 1:T_pre[i], k]);
      X_y_c[i, , k] = X_y[i, , k] - X_y_mean[i, k];
    }
  }
  // centered covariates for donors
  array[J] matrix[T, K] X_z_c;
  matrix[J, K] X_z_mean;
  {
    int min_T_pre = min(T_pre);
    for (j in 1:J) {
      for (k in 1:K) {
        X_z_mean[j, k] = mean(X_z[j, 1:min_T_pre, k]);
        X_z_c[j, , k] = X_z[j, , k] - X_z_mean[j, k];
      }
    }
  }
