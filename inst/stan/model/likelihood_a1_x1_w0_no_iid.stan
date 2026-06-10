  if (likelihood) {
    for (i in 1:N) {
      int Ti = T_pre[i];
      head(y[i], Ti) ~ normal(alpha[i] + X_y[i, 1:Ti] * beta, sigma[i]);
    }
  }
