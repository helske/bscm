  // center covariates using pre-treatment period
  // working intercept a_y absobs the shifts
  
  // centered covariates for treated
  matrix[T_pre, K] X_y_c = X_y[1:T_pre];
  row_vector[K] X_y_mean;
  for (k in 1:K) {
    X_y_mean[k] = mean(X_y_c[, k]);
    X_y_c[, k] -= X_y_mean[k];
  }
  // centered covariates for donors
  array[J] matrix[T_pre, K] X_z_c;
  matrix[J, K] X_z_mean;
  for (j in 1:J) {
    for (k in 1:K) {
      X_z_mean[j, k] = mean(X_z[j, 1:T_pre, k]);
      X_z_c[j, , k] = X_z[j, 1:T_pre, k] - X_z_mean[j, k];
    }
  }
