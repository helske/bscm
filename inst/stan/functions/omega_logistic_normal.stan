/**
   * Constrain sum-to-zero vector
   *
   * CRAN version of Stan does not support sum_to_zero_vector type
   * This is a workaround by Brian Ward to mimic this in older versions: 
   * https://discourse.mc-stan.org/t/timeline-for-rstan-update-on-cran/41038/2
   *
   * @param y unconstrained zero-sum parameters
   * @return vector z, the vector whose elements sum to zero
   */
  vector zero_sum_constrain(vector y) {
    int N = num_elements(y);
    vector[N + 1] z = zeros_vector(N + 1);
    real sum_w = 0;
    for (ii in 1:N) {
      int i = N - ii + 1; 
      real n = i;
      real w = y[i] * inv_sqrt(n * (n + 1));
      sum_w += w;
      z[i] += sum_w;     
      z[i + 1] -= w * n;    
    }
    return z;
  }
