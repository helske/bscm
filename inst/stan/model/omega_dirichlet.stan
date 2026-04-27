  for (i in 1:N) {
    omega[i] ~ dirichlet(kappavec);
  }
