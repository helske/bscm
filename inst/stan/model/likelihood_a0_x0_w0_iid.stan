  if (likelihood) {
    for (i in 1:N) {
      int Ti = T_pre[i];
      head(y[i], Ti) ~ normal(Z[1:Ti] * omega[i], sigma[i]);
    }
  }
