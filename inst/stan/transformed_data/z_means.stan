  array[N] row_vector[J] Z_mean;
  for (i in 1:N) {
    for (j in 1:J) {
      Z_mean[i, j] = mean(Z[1:T_pre[i], j]);
    }
  }
