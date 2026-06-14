    matrix[T, J] Z_term = cZ;
    for (j in 1:J) {
      Z_term[, j] -= cX_z[j] * beta;
    }
