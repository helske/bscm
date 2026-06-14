  if (likelihood) {
    for (i in 1:N) {
      int Ti = T_pre[i];
      head(y[i], Ti) ~ normal(a[i], sigma[i]);
    }
  }
