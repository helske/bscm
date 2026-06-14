  if (likelihood) {
    vector[T] mu;
    for (i in 1:N) {
      int Ti = T_pre[i];
      mu = a[i] + cZ * omega[i];
      mu[2:Ti] += rho[i] * head(y[i] - mu, Ti - 1);
      head(y[i], Ti) ~ normal(head(mu, Ti), sigma[i]);
    }
  }
