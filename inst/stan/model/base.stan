  for (i in 1:N) {
    eta[i] ~ std_normal();    
  }
  sigma ~ exponential(pr_rate_sigma);
