      y_mean[, i] = mu;
      y_mean[2:Ti, i] += rho[i] * head(y[i] - mu, Ti - 1);
      
      real r = y[i, Ti] - mu[Ti];
      for (t in (Ti + 1):T) {
        r *= rho[i];
        y_mean[t, i] = mu[t] + r;
      }
      if (cv == 0) {
        y_rep[1, i] = normal_rng(mu[1], sigma[i]);
        for (t in 2:T) {
          mu[t] += rho[i] * (y_rep[t - 1, i] - mu[t - 1]);
          y_rep[t, i] = normal_rng(mu[t], sigma[i]);
        }
      }
