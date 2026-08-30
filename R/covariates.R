#' Predictor arrays of the treated and donor units
#'
#' Returns the model matrix of the predictors as two arrays, one for the
#' treated units (`X_y`) and one for the donors (`X_z`), both with dimensions
#' unit, time, and predictor.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @noRd
get_Xs <- function(x) {
  X <- stats::model.matrix(stats::formula(x), data = x$data)
  if (has_intercept(x)) {
    X <- X[, -1L, drop = FALSE]
  }
  K <- ncol(X)
  T_total <- get_T_total(x)
  N <- get_N(x)
  J <- length(get_donors(x))
  n_units <- J + N
  X <- simplify2array(
    lapply(seq_len(K), \(k) matrix(X[, k], n_units, T_total, TRUE))
  )
  treated_idx <- which(unique(x$data[[get_unit(x)]]) %in% get_treated(x))
  X_z <- X[-treated_idx, , , drop = FALSE]
  X_y <- X[treated_idx, , , drop = FALSE]
  list(X_y = X_y, X_z = X_z)
}

#' Predictor differences between a treated unit and its synthetic control
#'
#' Returns the differences between the predictor values of the treated unit
#' `unit` and the corresponding donor weighted predictor values, as a time by
#' predictor `rvar` matrix. These differences are the covariate imbalances
#' before averaging over the predictors, and the multipliers of the regression
#' coefficients in the covariate adjustments.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param unit \[`integer(1)`]\cr Index of the treated unit.
#' @param X Output of `get_Xs()`.
#' @noRd
covariate_delta <- function(x, unit, X) {
  omega <- rvars_of(x, "omega")[, unit, drop = TRUE]
  do.call(
    cbind,
    lapply(
      seq_len(dim(X$X_y)[3L]),
      \(k) X$X_y[unit, , k] - c(omega %*% X$X_z[,, k])
    )
  )
}
