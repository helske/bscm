// Intercept-only model, no donors
data {
  int<lower = 1> J;                       // number of donor series
  int<lower = 1> N;                       // number of treated series
  int<lower = 1> T;                       // number of time points
  array[N] int<lower=1> T_pre;            // number of pre-treatment time points
  array[N] vector[T] y;                   // treated
  vector<lower = 0>[N] pr_rate_sigma;     // sigma ~ exponential()
  vector[N] pr_mean_intercept;                 // prior mean of intercept
  vector<lower = 0>[N] pr_sd_intercept;        // prior SD of intercept
}

parameters {
  vector<lower = 0>[N] sigma; // SD of the error term
  vector[N] a; // working intercept
}
transformed parameters {
  vector[N] alpha = a;
}
model {
  sigma ~ exponential(pr_rate_sigma);
  a ~ normal(pr_mean_intercept, pr_sd_intercept);
  for (i in 1:N) {
    int Ti = T_pre[i];
    head(y[i], Ti) ~ normal(a[i], sigma[i]);
  }
}
generated quantities {
  matrix[T, N] y_mean;
  for (i in 1:N) {
    y_mean[, i] = rep_vector(alpha[i], T);
  }
  matrix[T, N] y_rep;
  for (i in 1:N) {
    y_rep[, i] = to_vector(normal_rng(y_mean[, i], sigma[i]));
  }
}
