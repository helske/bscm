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
  array[N] vector[T] y;            // treated
  matrix[T, J] Z;                  // donors
  array[N] matrix[T, K] X_y;       // covariates for treated
  array[J] matrix[T, K] X_z;       // covariates for donors
  array[L] int<lower = 1, upper = K> tv_idx; // indices of gamma coefs
  int<lower = 0, upper = 1> sample_y_rep; // if 1, samples y_rep
}
transformed data {
  array[K] matrix[T, J] X_z_;
  if (K > 0) {
    for (j in 1:J) {
      for (k in 1:K) {
        X_z_[k, , j] = X_z[j, , k];
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
  matrix[J, N] omega;
  vector[N] sigma;                   // SD of the error term
  vector[N * icpt] alpha;                       // intercept
  vector[K] beta;                               // regression coefficients
  matrix[T, L] gamma;                           // time-varying coefficients
  vector[N * ar1] rho;   // autoregressive parameters
}
generated quantities {
  matrix[T, N] y_mean = Z * omega;
  matrix[T * sample_y_rep, N * sample_y_rep] y_rep;
  if (K > 0) {
    for (k in 1:K) {
      y_mean -= beta[k] * X_z_[k] * omega;
    }
    if (L > 0) {
      for (l in 1:L) {
        y_mean -= diag_pre_multiply(gamma[, l], W_z[l] * omega);
      }
    }
  }
  for (i in 1:N) {
    if (icpt) y_mean[, i] += alpha[i];
    if (K > 0) {
      y_mean[, i] += X_y[i] * beta;
      if (L > 0) y_mean[, i] += rows_dot_product(W_y[i], gamma);
    }
    if (ar1) {
      y_mean[2:T, i] += rho[i] * (y[i, 1:(T - 1)] - y_mean[1:(T - 1), i]);
    }
  }
  if (sample_y_rep) {
    for (i in 1:N) {
      y_rep[, i] = to_vector(normal_rng(y_mean[, i], sigma[i]));
    }
  }
}
