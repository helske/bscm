  //not supported on the current CRAN versions (as in April 2026)
  //array[N] sum_to_zero_vector[J] eta;
  array[N * dirichlet_omega] simplex[J] omega_;
  array[N * !dirichlet_omega] vector[J - 1] eta;
