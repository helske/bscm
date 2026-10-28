// marginal log-likelihood log p(y | mu) where y_i ~ N(mu_i, sigma^2) iid,
// i = 1, ..., n, with sigma integrated out, as a function of the residual sum
// of squares rss = ||y - mu||^2. dist is the code of the prior of sigma
// as in prior_lpdf():
//   5: half-normal
//   9: scaled inverse chi
//   10: square-root inverse gamma
real marginal_log_likelihood_lpdf(real rss, data int n, data int dist, data vector pars) {
  if (dist == 5) {
    return marginalized_normal_half_normal_lupdf(rss | n, pars[1]);
  } else if (dist == 9) {
    return marginalized_normal_scaled_inv_chi_lupdf(rss | n, pars[1], pars[2]);
  } else if (dist == 10) {
    return marginalized_normal_sqrt_inv_gamma_lupdf(rss | n, pars[1], pars[2]);
  }
  reject("Unknown prior distribution code: ", dist);
}

// Draw sigma from its conditional posterior given the residual sum of squares
// rss of n observations y_i ~ N(mu_i, sigma^2) iid, where dist is the code of
// the prior of sigma as above. With n = 0, sigma is drawn from its prior.
real marginalized_sigma_rng(real rss, int n, data int dist, data vector pars) {
  if (dist == 5) {
      if (n == 0) {
      return abs(normal_rng(0, pars[1]));
    }
    // with prior sigma ~ half-normal(s), the conditional posterior density of
    // sigma^2 is proportional to
    //   (sigma^2)^(-(n + 1) / 2) * exp(-(rss / sigma^2 + sigma^2 / s^2) / 2),
    // a generalized inverse Gaussian distribution GIG(lambda, psi, chi) with 
    // lambda = - (n - 1) / 2, psi = 1 / pars[1]^2, chi = rss. 
    return sqrt(gig_rng(-0.5 * (n - 1), inv_square(pars[1]), rss));
  } else if (dist == 9) {
    // sigma^2 ~ Inv-Gamma((v + n) / 2, (s^2 + rss) / 2)
    return sqrt(inv_gamma_rng(0.5 * (pars[1] + n), 0.5 * (square(pars[2]) + rss)));
  } else if (dist == 10) {
    // sigma^2 ~ Inv-Gamma(a + n / 2, b + rss / 2)
    return sqrt(inv_gamma_rng(pars[1] + 0.5 * n, pars[2] + 0.5 * rss));
  }
  reject("Unknown prior distribution code: ", dist);
}

// sigma ~ sqrt_inv_gamma(a, b)
real marginalized_normal_sqrt_inv_gamma_lpdf(real rss, int n, data real a, data real b) {
  return lgamma(a + 0.5 * n) - lgamma(a) - 0.5 * n * log(2 * pi() * b)
  - (a + 0.5 * n) * log1p(0.5 * rss / b);
}
// sigma ~ scaled-inverse-chi(v, s) <=> sigma ~ sqrt_inv_gamma(v/2, s^2/2)
real marginalized_normal_scaled_inv_chi_lpdf(real rss, int n, data real v, data real s) {
  real a = 0.5 * v;
  real b = 0.5 * square(s);

  return lgamma(a + 0.5 * n) - lgamma(a) - 0.5 * n * log(2 * pi() * b)
  - (a + 0.5 * n) * log1p(0.5 * rss / b);
}

// sigma ~ half-normal(s)
//   p(y) \propto r^{-(n-1)/2} K_{(n-1)/2}(r / s), where r = ||y - mu||
// Note: the density is infinite at r = 0).
real marginalized_normal_half_normal_lpdf(real rss, data int n, data real s) {
  return log2() - 0.5 * (n + 1) * log(2 * pi() * s) -
  0.25 * (n - 1) * log(rss) + log_bessel_K(n, sqrt(rss) / s);
}

// Generalized inverse Gaussian distribution GIG(lambda, psi, chi) with density
//   f(x | lambda, psi, chi) \propto x^(lambda - 1) exp(-(chi / x + psi x) / 2),
// for x > 0 and psi, chi > 0, sampled following Hörmann, W. and Leydold, J.
// (2013). Generating generalized inverse Gaussian random variates. Statistics
// and Computing 24, 547-557.
// With alpha = sqrt(psi / chi) and beta = sqrt(psi * chi), X = Y / alpha, where
// Y has the quasi-density
//   g(y | lambda, beta) = y^(lambda - 1) * exp(-beta / 2 (y + 1 / y)).
// For lambda < 0, 1 / X ~ GIG(-lambda, chi, psi), so that X = 1 / (alpha * Y),
// where Y ~ g(y | -lambda, beta). Y is drawn with the Ratio-of-Uniforms
// method with mode shift and minimal bounding rectangle (Algorithm 3), for
// which the rejection constant is at most 2 for |lambda| >= 1. For numerical
// stability, g is normalized by its value at the mode m, so that v+ = 1.
real gig_rng(real lambda, real psi, real chi) {
  real alpha = sqrt(psi / chi);
  real beta = sqrt(psi * chi);
  real lambda_abs = abs(lambda);
  real m = (sqrt(square(lambda_abs - 1) + square(beta)) + lambda_abs - 1) / beta;
  real log_gm = gig_log_quasi_density(m, lambda_abs, beta);
  // roots x- in (0, m) and x+ in (m, Inf) of the cubic
  //   x^3 + a x^2 + b x + c = 0
  // using Cardano's formula
  real a = -2 * (lambda_abs + 1) / beta - m;
  real b = 2 * (lambda_abs - 1) * m / beta - 1;
  real c = m;
  real p = b - square(a) / 3;
  real q = 2 * pow(a, 3) / 27 - a * b / 3 + c;
  real phi = acos(-q / 2 * sqrt(-27 / pow(p, 3)));
  real x_minus = sqrt(-4 * p / 3) * cos(phi / 3 + 4 * pi() / 3) - a / 3;
  real x_plus = sqrt(-4 * p / 3) * cos(phi / 3) - a / 3;
  real u_minus = (x_minus - m) *
    exp(0.5 * (gig_log_quasi_density(x_minus, lambda_abs, beta) - log_gm));
  real u_plus = (x_plus - m) *
    exp(0.5 * (gig_log_quasi_density(x_plus, lambda_abs, beta) - log_gm));
  // the expected number of iterations is at most 2 for |lambda| >= 1, the
  // upper limit guards against an infinite loop due to numerical problems
  int i = 0;
  while (i < 1000) {
    real u = uniform_rng(u_minus, u_plus);
    real v = uniform_rng(0, 1);
    real x = u / v + m;
    if (x > 0 && 2 * log(v) <= gig_log_quasi_density(x, lambda_abs, beta) - log_gm) {
      return lambda < 0 ? 1 / (alpha * x) : x / alpha;
    }
    i += 1;
  }
  reject("No GIG draw accepted in 1000 iterations with lambda = ",
         lambda, ", psi = ", psi, ", chi = ", chi);
}
real gig_log_quasi_density(real x, real lambda, real beta) {
  return (lambda - 1) * log(x) - 0.5 * beta * (x + 1 / x);
}

// Bessel function stuff for halfnormal prior marginalization

// log K_nu(x) for nu in {0, 1}, stable for large x.
// Below x = 500 uses Stan's built-in Bessel; above it, the large-x
// asymptotic expansion (4 terms, relative error < 1e-15 there).
real log_bessel_K_01(data int nu, real x) {
  if (x < 500) {
    return log(modified_bessel_second_kind(nu, x));
  }
  real mu4 = 4 * nu^2;
  real t = 1;
  real term = 1;
  for (k in 1:4) {
    term *= (mu4 - (2 * k - 1)^2) / (k * 8 * x);
    t += term;
  }
  return 0.5 * log(pi() / (2 * x)) - x + log(t);
}
// log K_{(n-1)/2}(x), computed on the log scale by upward recurrence
//   K_{v+1}(x) = K_{v-1}(x) + (2 v / x) K_v(x)
// starting from K_0, K_1 (n odd) or K_{1/2}, K_{3/2} (n even).
real log_bessel_K(data int n, real x) {
  real nu0;
  int m;
  real a;  // log K_{nu0 + k - 1}
  real b;  // log K_{nu0 + k}
  if (n % 2 == 1) {
    nu0 = 0;
    m = (n - 1) %/% 2;
    a = log_bessel_K_01(0, x);
    b = log_bessel_K_01(1, x);
  } else {
    nu0 = 0.5;
    m = n %/% 2 - 1;
    a = 0.5 * log(pi() / (2 * x)) - x;  // K_{1/2}, exact
    b = a + log1p(1 / x);               // K_{3/2}, exact
  }
  if (m == 0) {
    return a;
  }
  for (k in 1:(m - 1)) {
    real b_new = log_sum_exp(a, log(2 * (nu0 + k) / x) + b);
    a = b;
    b = b_new;
  }
  return b;
}

