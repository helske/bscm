functions {
  /*
  * Constrain sum-to-zero vector
  *
  * CRAN version of Stan does not support sum_to_zero_vector type
  * This is a workaround by Brian Ward to mimic this in older versions:
  * https://discourse.mc-stan.org/t/timeline-for-rstan-update-on-cran/41038/2
  *
  * @param y unconstrained zero-sum parameters
  * @return vector z, the vector whose elements sum to zero
  */
  vector zero_sum_constrain(vector y) {
    int N = num_elements(y);
    vector[N + 1] z = zeros_vector(N + 1);
    real sum_w = 0;
    for (ii in 1 : N) {
      int i = N - ii + 1;
      real n = i;
      real w = y[i] * inv_sqrt(n * (n + 1));
      sum_w += w;
      z[i] += sum_w;
      z[i + 1] -= w * n;
    }
    return z;
  }
}
data {
  int<lower=0, upper=1> use_alpha;
  int<lower=0, upper=1> use_beta;
  int<lower=0, upper=1> use_gamma;
  int<lower=0, upper=1> use_ar1;
  int<lower=2> J; // number of donor series
  int<lower=1> N; // number of treated series
  int<lower=2> T; // number of time points
  array[N] int<lower=1, upper=T> T_pre; // number of pre-treatment time points
  int<lower=use_beta> K; // number of covariates
  int<lower=use_gamma, upper=K> L; // number of time-varying coefficients
  int<lower=1> D; // number of spline coefficients + 1
  array[N] vector[T] y; // treated
  matrix[T, J] Z; // donors
  array[use_beta ? N : 0] matrix[T, K] X_y; // covariates for treated
  array[use_beta ? J : 0] matrix[T, K] X_z; // covariates for donors
  array[L] int<lower=1, upper=K> tv_idx; // indices of gamma coefs
  matrix[T, D - 1] A; // spline projection matrix
  int<lower=0, upper=1> likelihood; // if 0, sample from prior only
  int<lower=0, upper=8> pr_dist_sigma;
  int<lower=0, upper=8> pr_dist_omega;
  int<lower=0, upper=8> pr_dist_intercept;
  int<lower=0, upper=8> pr_dist_beta;
  int<lower=0, upper=8> pr_dist_kappa;
  int<lower=0, upper=8> pr_dist_rho;
  matrix[N, 2] pr_pars_sigma;
  vector[N] pr_pars_omega;
  matrix[use_alpha ? N : 0, 3] pr_pars_intercept;
  matrix[K, 3] pr_pars_beta;
  matrix[L, 2] pr_pars_kappa;
  matrix[use_ar1 ? N : 0, 2] pr_pars_rho;
}
transformed data {
  int T0 = min(T_pre);
  
  // vector of concentration parameters for Dirichlet prior of omega
  // or standard deviations for logistic normal
  int<lower=0, upper=1> dirichlet_omega = pr_dist_omega == 7;
  matrix[(J - 1) * dirichlet_omega + 1, N] pr_omega;
  // scaling for logistic normal prior of omega
  if (dirichlet_omega) {
    for (i in 1 : N) {
      pr_omega[ : , i] = rep_vector(pr_pars_omega[i], J);
    }
  } else {
    // adjust the marginal standard deviation of elements of sum-to-zero eta
    // eta ~ N(0, scale), sum(eta) = 0 => marginal SD(eta_i) = 1
    pr_omega[1,  : ] = pr_pars_omega' * sqrt(J * inv(J - 1));
  }
  
  // centered donors
  row_vector[J] Z_mean = rep_row_vector(0, J);
  matrix[T, J] cZ = Z;
  if (use_alpha) {
    for (j in 1 : J) {
      Z_mean[j] = mean(Z[1 : T0, j]);
      cZ[ : , j] -= Z_mean[j];
    }
  }
  
  matrix[N, K] X_y_mean;
  matrix[J, K] X_z_mean;
  array[N] matrix[T, K] cX_y;
  array[K] matrix[T, J] cX_z;
  if (use_beta) {
    for (i in 1 : N) {
      for (k in 1 : K) {
        X_y_mean[i, k] = mean(X_y[i, 1 : T0, k]);
        cX_y[i,  : , k] = X_y[i,  : , k] - X_y_mean[i, k];
      }
    }
    for (j in 1 : J) {
      for (k in 1 : K) {
        X_z_mean[j, k] = mean(X_z[j, 1 : T0, k]);
        cX_z[k,  : , j] = X_z[j,  : , k] - X_z_mean[j, k];
      }
    }
  }
  array[N] matrix[T, L] W_y;
  array[L] matrix[T, J] W_z;
  if (use_gamma) {
    for (i in 1 : N) {
      W_y[i] = X_y[i,  : , tv_idx];
    }
    for (l in 1 : L) {
      for (j in 1 : J) {
        W_z[l,  : , j] = X_z[j,  : , tv_idx[l]];
      }
    }
  }
}
parameters {
  //not supported on the current CRAN versions (as in April 2026)
  //array[N] sum_to_zero_vector[J] eta;
  array[N * dirichlet_omega] simplex[J] omega_; // omega ~ dirichlet
  array[N * !dirichlet_omega] vector[J - 1] eta; // omega ~ logistic normal
  vector<lower=0>[N] sigma; // SD of the error term
  vector[use_alpha ? N : 0] a; // working intercept
  vector[use_beta ? K : 0] beta; // regression coefficients
  matrix[D - 1, use_gamma ? L : 0] xi; // spline coefficients
  vector<lower=0>[use_gamma ? L : 0] kappa; // SD of RW1 prior xi
  vector<lower=-1, upper=1>[use_ar1 ? N : 0] rho; // autoregressive parameters
}
transformed parameters {
  matrix[J, N] omega;
  if (dirichlet_omega) {
    for (i in 1 : N) {
      omega[ : , i] = omega_[i];
    }
  } else {
    for (i in 1 : N) {
      omega[ : , i] = softmax(pr_omega[1, i] * zero_sum_constrain(eta[i]));
    }
  }
  // time-varying regression coefficients
  matrix[use_gamma ? T : 0, L] gamma;
  for (l in 1 : L) {
    gamma[ : , l] = A * (kappa[l] * xi[ : , l]);
  }
}
model {
  // prior for sigma
  if (pr_dist_sigma == 3) {
    sigma ~ exponential(pr_pars_sigma[ : , 1]);
  } else if (pr_dist_sigma == 4) {
    sigma ~ gamma(pr_pars_sigma[ : , 1], pr_pars_sigma[ : , 2]);
  } else {
    sigma ~ normal(0, pr_pars_sigma[ : , 1]);
  }
  // prior for omega
  if (dirichlet_omega) {
    for (i in 1 : N) {
      omega_[i] ~ dirichlet(pr_omega[ : , i]);
    }
  } else {
    for (i in 1 : N) {
      eta[i] ~ std_normal();
    }
  }
  // prior for alpha
  if (pr_dist_intercept == 1) {
    a ~ normal(pr_pars_intercept[ : , 1], pr_pars_intercept[ : , 2]);
  } else if (pr_dist_intercept == 2) {
    a ~ student_t(pr_pars_intercept[ : , 1], pr_pars_intercept[ : , 2],
                  pr_pars_intercept[ : , 3]);
  }
  // prior for beta
  if (pr_dist_beta == 1) {
    beta ~ normal(pr_pars_beta[ : , 1], pr_pars_beta[ : , 2]);
  } else if (pr_dist_beta == 2) {
    beta ~ student_t(pr_pars_beta[ : , 1], pr_pars_beta[ : , 2],
                     pr_pars_beta[ : , 3]);
  }
  // prior for gamma
  if (pr_dist_kappa == 3) {
    kappa ~ exponential(pr_pars_kappa[ : , 1]);
  } else if (pr_dist_kappa == 4) {
    kappa ~ gamma(pr_pars_kappa[ : , 1], pr_pars_kappa[ : , 2]);
  } else if (pr_dist_kappa == 6) {
    kappa ~ normal(0, pr_pars_kappa[ : , 1]);
  }
  to_vector(xi) ~ std_normal();
  // prior for rho
  if (pr_dist_rho == 5) {
    0.5 * (1 + rho) ~ beta(pr_pars_rho[ : , 1], pr_pars_rho[ : , 2]);
  }
  if (likelihood) {
    matrix[T, N] Z_term = cZ * omega;
    if (K > 0) {
      for (k in 1 : K) {
        Z_term -= beta[k] * cX_z[k] * omega;
      }
      if (L > 0) {
        for (l in 1 : L) {
          Z_term -= diag_pre_multiply(gamma[ : , l], W_z[l] * omega);
        }
      }
    }
    for (i in 1 : N) {
      int Ti = T_pre[i];
      vector[Ti] mu = Z_term[1 : Ti, i];
      if (use_alpha) {
        mu += a[i];
      }
      if (use_beta) {
        mu += cX_y[i, 1 : Ti,  : ] * beta;
      }
      if (use_gamma) {
        mu += rows_dot_product(W_y[i, 1 : Ti,  : ], gamma[1 : Ti,  : ]);
      }
      if (use_ar1) {
        mu[2 : Ti] += rho[i] * (y[i, 1 : (Ti - 1)] - mu[1 : (Ti - 1)]);
      }
      y[i, 1 : Ti] ~ normal(mu, sigma[i]);
    }
  }
}
generated quantities {
  vector[use_alpha ? N : 0] alpha = rep_vector(0, N * use_alpha);
  if (use_alpha) {
    alpha = a;
    for (i in 1 : N) {
      alpha[i] -= Z_mean * omega[ : , i];
    }
    if (K > 0) {
      alpha -= X_y_mean * beta;
      vector[J] xzb = X_z_mean * beta;
      for (i in 1 : N) {
        alpha[i] += dot_product(xzb, omega[ : , i]);
      }
    }
  }
}
