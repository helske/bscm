real sqrt_inv_gamma_lpdf(vector y, vector alpha, vector beta) {
  return inv_gamma_lupdf(square(y) | alpha, beta) + sum(log(2 * y));
}

real scaled_inv_chi_lpdf(vector y, vector nu, vector sigma) {
  return scaled_inv_chi_square_lupdf(square(y) | nu, sigma .* inv_sqrt(nu)) + 
  sum(log(2 * y));
}
/**
* Prior log-densities for scalar parameters (vectorized)
* Codes (alphabetical order, matching prior_to_stan() in R):
*   0: none
*   1: beta
*   2: cauchy
*   3: exponential
*   4: gamma
*   5: half-normal
*   6: inverse gamma
*   7: laplace (double_exponential)
*   8: normal
*   9: scaled inverse chi
*   10: square-root inverse gamma
*   11: t
*/
real prior_lpdf(vector x, data int dist, data matrix pars) {
  if (dist == 0) {
    return 0;
  } else if (dist == 1) {
    return beta_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 2) {
    return cauchy_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 3) {
    return exponential_lupdf(x | pars[, 1]);
  } else if (dist == 4) {
    return gamma_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 5) {
    return normal_lupdf(x | 0, pars[, 1]);
  } else if (dist == 6) {
    return inv_gamma_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 7) {
    return double_exponential_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 8) {
    return normal_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 9) {
    return scaled_inv_chi_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 10) {
    return sqrt_inv_gamma_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 11) {
    return student_t_lupdf(x | pars[, 1], pars[, 2], pars[, 3]);
  }
  reject("Unknown prior distribution code: ", dist);
}
// Overload for dist == 0
real prior_lpdf(vector x, int dist) {
  return prior_lupdf(x | 0, rep_matrix(0, 0, 0));
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
