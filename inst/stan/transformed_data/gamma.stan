  row_vector[L] L_zeros = rep_row_vector(0, L);
  array[N] matrix[T, L] W_y;
  for (i in 1:N) {
    W_y[i] = X_y[i, , tv_idx];
  }
  array[J] matrix[T, L] W_z;
  for (j in 1:J) {
    W_z[j] = X_z[j, , tv_idx];
  }
