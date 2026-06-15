  if (dirichlet_omega) {
    for (i in 1:N) {
      omega_[i] ~ dirichlet(kappavec);
    }
  } else {
    for (i in 1:N) {
      eta[i] ~ std_normal();
    }
  }
