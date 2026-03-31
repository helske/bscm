  for (i in 1:N) {
    omega[i] ~ dirichlet(kappavec);
  }
  sigma ~ exponential(pr_rate_sigma);
