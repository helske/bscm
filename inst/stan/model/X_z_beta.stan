    matrix[T, J] Z_term = Z;
    for (j in 1:J) {
      Z_term[, j] -= X_z[j] * beta;
    }
