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
    for (ii in 1:N) {
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
  int<lower = 1> J;                // number of donor series
  int<lower = 1> N;                // number of treated series
  int<lower = 1> T;                // number of time points
  array[N] int<lower=1> T_pre;     // number of pre-treatment time points
  int<lower = 0> K;                // number of covariates
  int<lower = 0, upper = K> L;     // number of time-varying coefficients
  int<lower = 1> D;                // number of spline coefficients + 1
  int<lower = 0, upper = 1> icpt;  // 0 = no intercept, 1 = intercept
  int<lower = 0, upper = 1> ar1;   // 0 = iid errors, 1 = AR(1) errors
  int<lower = 0, upper = 1> dirichlet_omega; // dirichlet or logistic normal
  real<lower = 0> kappa;           // hyperparameter of prior on omega
  array[N] vector[T] y;            // treated
  matrix[T, J] Z;                  // donors
  array[N] matrix[T, K] X_y;       // covariates for treated
  array[J] matrix[T, K] X_z;       // covariates for donors
  array[L] int<lower = 1, upper = K> tv_idx;    // indices of gamma coefs
  matrix[T, D - 1] A;                           // spline projection matrix
  vector<lower = 0>[N] pr_rate_sigma;           // sigma ~ exponential()
  vector[N * icpt] pr_mean_intercept;           // prior mean of intercept
  vector<lower = 0>[N * icpt] pr_sd_intercept;  // prior SD of intercept
  vector[K] pr_mean_beta;                       // prior mean of beta
  vector<lower = 0>[K] pr_sd_beta;              // prior SD of beta
  vector<lower = 0>[L] pr_sd_sigma_gamma;
  vector<lower = 0>[N * ar1] pr_shape1_rho;     // 0.5 * (1 + rho) ~ Beta()
  vector<lower = 0>[N * ar1] pr_shape2_rho;
  int<lower = 0, upper = 1> likelihood;         // if 0, sample from prior only
}
transformed data {
  int T0 = min(T_pre);
  // vector of concentration parameters for Dirichlet prior of omega
  // standard deviations for logistic normal
  
  vector[(J - 1) * dirichlet_omega + 1] kappavec;
  // scaling for logistic normal prior of omega
  if (dirichlet_omega) {
    kappavec = rep_vector(kappa, J);
  } else {
    // adjust the marginal standard deviation of elements of sum-to-zero eta
    // eta ~ N(0, scale), sum(eta) = 0 => marginal SD(eta_i) = 1
    kappavec[1] = kappa * sqrt(J * inv(J - 1));
  }
  
  // centered donors
  row_vector[J] Z_mean = rep_row_vector(0, J);
  matrix[T, J] cZ = Z;
  if (icpt) {
    for (j in 1:J) {
      Z_mean[j] = mean(Z[1:T0, j]);
      cZ[, j] -= Z_mean[j];
    }
  }
  
  matrix[N, K] X_y_mean;
  matrix[J, K] X_z_mean;
  array[N] matrix[T, K] cX_y;
  array[K] matrix[T, J] cX_z;
  if (K > 0) {
    for (i in 1:N) {
      for (k in 1:K) {
        X_y_mean[i, k] = mean(X_y[i, 1:T0, k]);
        cX_y[i, , k] = X_y[i, , k] - X_y_mean[i, k];
      }
    }
    for (j in 1:J) {
      for (k in 1:K) {
        X_z_mean[j, k] = mean(X_z[j, 1:T0, k]);
        cX_z[k, , j] = X_z[j, , k] - X_z_mean[j, k];
      }
    }
  }
  array[N] matrix[T, L] W_y;
  array[L] matrix[T, J] W_z;
  if (L > 0) {
    for (i in 1:N) {
      W_y[i] = X_y[i, , tv_idx];
    }
    for (l in 1:L) {
      for (j in 1:J) {
        W_z[l, , j] = X_z[j, , tv_idx[l]];
      }
    }
  }
}
parameters {
  //not supported on the current CRAN versions (as in April 2026)
  //array[N] sum_to_zero_vector[J] eta;
  array[N * dirichlet_omega] simplex[J] omega_;
  array[N * !dirichlet_omega] vector[J - 1] eta;
  vector<lower = 0>[N] sigma;                   // SD of the error term
  vector[N * icpt] a;                           // working intercept
  vector[K] beta;                               // regression coefficients
  matrix[D - 1, L] xi;                          // spline coefficients
  vector<lower = 0>[L] sigma_gamma;             // SD of RW1 prior xi
  vector<lower = -1, upper = 1>[N * ar1] rho;   // autoregressive parameters
}
transformed parameters {
  matrix[J, N] omega;
  if (dirichlet_omega) {
    for (i in 1:N) {
      omega[, i] = omega_[i];
    }
  } else {
    for (i in 1:N) {
      omega[, i] = softmax(kappavec[1] * zero_sum_constrain(eta[i]));
    }
  }
  // time-varying regression coefficients
  matrix[T, L] gamma;
  for (l in 1:L) {
    gamma[, l] = A * (sigma_gamma[l] * xi[, l]);
  }
}
model {
  sigma ~ exponential(pr_rate_sigma);
  if (dirichlet_omega) {
    for (i in 1:N) {
      omega_[i] ~ dirichlet(kappavec);
    }
  } else {
    for (i in 1:N) {
      eta[i] ~ std_normal();
    }
  }
  a ~ normal(pr_mean_intercept, pr_sd_intercept);
  beta ~ normal(pr_mean_beta, pr_sd_beta);
  sigma_gamma ~ normal(0, pr_sd_sigma_gamma);
  to_vector(xi) ~ std_normal();
  0.5 * (1 + rho) ~ beta(pr_shape1_rho, pr_shape2_rho);
  if (likelihood) {
    matrix[T, N] Z_term = cZ * omega;
    if (K > 0) {
      for (k in 1:K) {
        Z_term -= beta[k] * cX_z[k] * omega;
      }
      if (L > 0) {
        for (l in 1:L) {
          Z_term -= diag_pre_multiply(gamma[, l], W_z[l] * omega);
        }
      }
    }
    for (i in 1:N) {
      int Ti = T_pre[i];
      vector[Ti] mu = Z_term[1:Ti, i];
      if (icpt) mu += a[i];
      if (K > 0) {
        mu += cX_y[i, 1:Ti, ] * beta;
        if (L > 0) mu += rows_dot_product(W_y[i, 1:Ti, ], gamma[1:Ti, ]);
      }
      if (ar1) mu[2:Ti] += rho[i] * (y[i, 1:(Ti - 1)] - mu[1:(Ti - 1)]);
      y[i, 1:Ti] ~ normal(mu, sigma[i]);
    }
  }
}
generated quantities {
  vector[N * icpt] alpha = rep_vector(0, N * icpt);
  if (icpt) {
    alpha = a;
    for (i in 1:N) {
      alpha[i] -= Z_mean * omega[, i];
    }
    if (K > 0) {
      alpha -= X_y_mean * beta;
      vector[J] xzb = X_z_mean * beta;
      for (i in 1:N) {
        alpha[i] += dot_product(xzb, omega[, i]);
      }
    }
  }
}
