  if (likelihood) {
    for (i in 1:N) {
      int Ti = T_pre[i];
      vector[Ti] mu = a[i] + cX_y[i, 1:Ti] * beta;
      mu[2:Ti] += rho[i] * (head(y[i] - mu, Ti - 1) - head(mu, Ti - 1));
      head(y[i], Ti) ~ normal(mu, sigma[i]);
    }
  }
