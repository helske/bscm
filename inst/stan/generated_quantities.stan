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
  int<lower=0, upper=1> sample_y_rep; // if 1, samples y_rep
}
transformed data {
  array[K] matrix[T, J] X_z_;
  if (use_beta) {
    for (j in 1 : J) {
      for (k in 1 : K) {
        X_z_[k,  : , j] = X_z[j,  : , k];
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
  matrix[J, N] omega; // donor weights
  vector[N] sigma; // SD of the error term
  vector[use_alpha ? N : 0] alpha; // intercept
  vector[use_beta ? K : 0] beta; // regression coefficients
  matrix[use_gamma ? T : 0, L] gamma; // time-varying regression coefficients
  vector[use_ar1 ? N : 0] rho; // autoregressive parameters
}
generated quantities {
  matrix[T, N] y_mean = Z * omega;
  matrix[T * sample_y_rep, N * sample_y_rep] y_rep;
  if (use_beta) {
    for (k in 1 : K) {
      y_mean -= beta[k] * X_z_[k] * omega;
    }
    if (use_gamma) {
      for (l in 1 : L) {
        y_mean -= diag_pre_multiply(gamma[ : , l], W_z[l] * omega);
      }
    }
  }
  for (i in 1 : N) {
    if (use_alpha) {
      y_mean[ : , i] += alpha[i];
    }
    if (use_beta) {
      y_mean[ : , i] += X_y[i] * beta;
    }
    if (use_gamma) {
      y_mean[ : , i] += rows_dot_product(W_y[i], gamma);
    }
    if (use_ar1) {
      int Ti = T_pre[i];
      y_mean[2 : Ti, i] += rho[i]
                           * (y[i, 1 : (Ti - 1)] - y_mean[1 : (Ti - 1), i]);
      real u = y[i, Ti] - y_mean[Ti, i];
      for (t in (Ti + 1) : T) {
        u *= rho[i];
        y_mean[t, i] += u;
      }
    }
  }
  if (sample_y_rep) {
    if (use_ar1) {
      for (i in 1 : N) {
        int Ti = T_pre[i];
        y_rep[1 : Ti, i] = to_vector(normal_rng(y_mean[1 : Ti, i], sigma[i]));
        real r = 0;
        for (t in (Ti + 1) : T) {
          r = normal_rng(rho[i] * r, sigma[i]);
          y_rep[t, i] = y_mean[t, i] + r;
        }
      }
    } else {
      for (i in 1 : N) {
        y_rep[ : , i] = to_vector(normal_rng(y_mean[ : , i], sigma[i]));
      }
    }
  }
}
