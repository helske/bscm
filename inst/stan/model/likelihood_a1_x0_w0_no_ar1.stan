  if (likelihood) {
    for (i in 1:N) {
      int Ti = T_pre[i];
      vector[Ti] mu = rep_vector(a[i], Ti);
      mu[2:Ti] += rho[i] * (head(y[i], Ti - 1) - a[i]);
      head(y[i], Ti) ~ normal(mu, sigma[i]);
    }
  }
