
// ILR transformation and symmetric Dirichlet(alpha) log-density using log_x
// dropping constant 0.5 * log(N)
vector dirichlet_ilr_simplex_lp(vector y, real alpha) {
  int N = rows(y) + 1;
  vector[N] z = zero_sum_constrain(y);
  real r = log_sum_exp(z);
  // jacobian from ILR transformation: const + sum(log_x), log_x = x - z
  // sum(log_x) = -N * r as sum(z) = 0
  // Dirichlet log-density:
  // (alpha - 1) * sum(log_x) + lgamma(N * alpha) - N * lgamma(alpha)
  target += - alpha * N * r + lgamma(N * alpha) - N * lgamma(alpha);
  return exp(z - r);
}

// with a tilt phi
vector dirichlet_ilr_simplex_lp(vector y, real alpha, vector phi) {
  int N = rows(y) + 1;
  vector[N] z = zero_sum_constrain(y);
  real r = log_sum_exp(z);
  target += - alpha * N * r + lgamma(N * alpha) - N * lgamma(alpha);
  return softmax(phi + z);  // = softmax(phi + log_x), r cancels
}

/* Create weights omega from unconstrained vector x and modify target accordingly
* Codes (alphabetical order, matching prior_to_stan() in R):
*  -1: dirichlet
*  -2: logistic laplace
*  -3: logistic normal
*/
vector build_omega_lp(vector x, real zeta, data int dist) {
  int J = rows(x) + 1;
  vector[J] omega;
  if (dist == -1) {
    omega = dirichlet_ilr_simplex_lp(x, zeta);
  } else {
    real scale = sqrt(J * inv(J - 1));
    if (dist == -2) {
      target += double_exponential_lpdf(x | 0, 1);
      scale *= inv_sqrt(2);
    } else if (dist == -3) {
      target += std_normal_lpdf(x);
    }
    omega = softmax(inv_sqrt(zeta) * scale * zero_sum_constrain(x));
  }
  return omega;
}
// As above but with a tilt phi
vector build_omega_lp(vector x, real zeta, data int dist, vector phi) {
  int J = rows(x) + 1;
  vector[J] omega;
  if (dist == -1) {
    omega = dirichlet_ilr_simplex_lp(x, zeta, phi);
  } else {
    real scale = sqrt(J * inv(J - 1));
    if (dist == -2) {
      target += double_exponential_lpdf(x | 0, 1);
      scale *= inv_sqrt(2);
    } else if (dist == -3) {
      target += std_normal_lpdf(x);
    }
    omega = softmax(phi + inv_sqrt(zeta) * scale * zero_sum_constrain(x));
  }
  return omega;
}
