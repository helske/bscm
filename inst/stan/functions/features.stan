// Reimplementations of some Stan features not available on version used by rstan
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
  vector[N + 1] z = rep_vector(0, N + 1);
  real sum_w = 0;
  for (ii in 1:N) {
    int i = N - ii + 1;
    real w = y[i] * inv_sqrt(i * (i + 1));
    sum_w += w;
    z[i] += sum_w;
    z[i + 1] -= w * i;
  }
  return z;
}

// ILR transformation and symmetric Dirichlet(alpha) log-density using log_x
// dropping constant 0.5 * log(N)
vector dirichlet_ilr_simplex_lp(vector y, real alpha) {
  int N = rows(y) + 1;
  vector[N] z = zero_sum_constrain(y);
  real r = log_sum_exp(z);
  vector[N] log_x = z - r;
  // jacobian from ILR transformation: const + sum(log_x)
  // sum(log_x) = -N * r as sum(z) = 0
  // Dirichlet log-density:
  // (alpha - 1) * sum(log_x) + lgamma(N * alpha) - N * lgamma(alpha)
  target += - alpha * N * r + lgamma(N * alpha) - N * lgamma(alpha);
  return exp(log_x);
}
