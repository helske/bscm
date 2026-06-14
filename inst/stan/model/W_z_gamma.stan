    for (j in 1:J) {
      Z_term[, j] -= rows_dot_product(W_z[j], gamma);
    }
