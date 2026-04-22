  sigma_gamma ~ gamma(2, pr_rate_sigma_gamma);
  to_vector(gamma_raw) ~ std_normal();
