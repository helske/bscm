  if (likelihood) {
    for (i in 1:N) {
      int Ti = T_pre[i];
      vector[Ti] mu = alpha[i] + X_y[i] * beta + 
      rows_dot_product(X_y[i, , tv_idx], gamma);
      head(y[i], Ti) ~ normal(head(mu, Ti), sigma[i]);
    }
  }
