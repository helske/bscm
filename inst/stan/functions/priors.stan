/**
* Prior log-densities for scalar parameters (vectorized)
* Codes:
*   0: none
*   1: normal
*   2: student_t
*   3: laplace (double_exponential)
*   4: exponential
*   5: gamma
*   6: half_normal
*   7: beta
*   8: cauchy
*/
real prior_lpdf(vector x, data int dist, data matrix pars) {
  if (dist == 0) {
    return 0;
  } else if (dist == 1) {
    return normal_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 2) {
    return student_t_lupdf(x | pars[, 1], pars[, 2], pars[, 3]);
  } else if (dist == 3) {
    return double_exponential_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 4) {
    return exponential_lupdf(x | pars[, 1]);
  } else if (dist == 5) {
    return gamma_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 6) {
    return normal_lupdf(x | 0, pars[, 1]);
  } else if (dist == 7) {
    return beta_lupdf(x | pars[, 1], pars[, 2]);
  } else if (dist == 8) {
    return cauchy_lupdf(x | pars[, 1], pars[, 2]);
  }
  reject("Unknown prior distribution code: ", dist);
}
// Overload for dist == 0
real prior_lpdf(vector x, int dist) {
  return prior_lupdf(x | 0, rep_matrix(0, 0, 0));
}



/* Create weights omega from unconstrained vector x and modify target accordingly
*  1: dirichlet
*  2: logistic normal
*  3: logistic laplace
*  4: logistic t
*/
vector build_omega_lp(vector x, real zeta, data int dist) {
  int J = rows(x) + 1;
  vector[J] omega;
  if (dist == -1) {
    omega = dirichlet_ilr_simplex_lp(x, zeta);
  } else {
    real scale = sqrt(J * inv(J - 1));
    if (dist == -2) {
      target += std_normal_lpdf(x);
    } else if (dist == -3) {
      target += double_exponential_lpdf(x | 0, 1);
      scale *= inv_sqrt(2);
    } else {
      target += student_t_lpdf(x | 3, 0, 1);
      scale *= inv_sqrt(3);
    }
    omega = softmax(zeta * scale * zero_sum_constrain(x));
  }
  return omega;
}
