  // adjust the marginal standard deviation of elements of sum-to-zero eta
  // eta ~ N(0, scale), sum(eta) = 0 => marginal SD(eta_i) = 1
  real<lower = 0> scale = sqrt(J * inv(J - 1));
