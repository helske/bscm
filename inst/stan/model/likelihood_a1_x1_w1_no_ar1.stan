  if (likelihood) {
    for (i in 1:N) {
      int Ti = T_pre[i];
      vector[Ti] mu = alpha[i] + X_y[i, 1:Ti] * beta +
      rows_dot_product(X_y[i, 1:Ti, tv_idx], gamma[1:Ti]);
      mu[2:Ti] += rho[i] * (head(y[i], Ti - 1) - head(mu, Ti - 1));
      head(y[i], Ti) ~ normal(mu, sigma[i]);
    }
  }
