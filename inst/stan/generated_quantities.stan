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
  array[N] matrix[T, L] W_y;
  array[J] matrix[T, L] W_z;
  if (L > 0) {
    for (i in 1:N) {
      W_y[i] = X_y[i, , tv_idx];
    }
    for (j in 1:J) {
      W_z[j] = X_z[j, , tv_idx];
    }
  }
}
parameters {
  array[N] vector[J] omega;
  vector<lower = 0>[N] sigma;                   // SD of the error term
  vector[N * icpt] alpha;                       // intercept
  vector[K] beta;                               // regression coefficients
  matrix[T, L] gamma;                           // time-varying coefficients
  vector<lower = -1, upper = 1>[N * ar1] rho;   // autoregressive parameters
}
generated quantities {
  matrix[T, N] y_mean;
  matrix[T * sample_y_rep, N * sample_y_rep] y_rep;
  {
    matrix[T, J] Z_term = Z;
    if (K > 0) {
      for (j in 1:J) {
        Z_term[, j] -= X_z[j] * beta;
      }
      if (L > 0) {
        for (j in 1:J) {
          Z_term[, j] -= rows_dot_product(W_z[j], gamma);
        }
      }
    }
    for (i in 1:N) {
      y_mean[, i] = Z_term * omega[i];
      if (icpt) y_mean[, i] += alpha[i];
      if (K > 0) {
        y_mean[, i] += X_y[i] * beta;
        if (L > 0) y_mean[, i] += rows_dot_product(W_y[i], gamma);
      }
      if (ar1) {
        y_mean[2:T, i] += rho[i] * (y[i, 1:(T - 1)] - y_mean[1:(T - 1), i]);
      }
      y_rep[, i] = to_vector(normal_rng(y_mean[, i], sigma[i]));
    }
  }
  if (sample_y_rep) {
    for (i in 1:N) {
      y_rep[, i] = to_vector(normal_rng(y_mean[, i], sigma[i]));
    }
  }
}
