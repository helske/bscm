#' Get input data to Stan from `bscmfit` object
#'
#' Reconstructs the data list passed to Stan in [bscm()]. Used by methods such
#' as [bscm::loo()] with `reloo = TRUE`.Requires the original data
#' (`save_data = TRUE`).
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A named list suitable for passing to [rstan::sampling()].
#' @export
get_standata <- function(x, ...) {
  UseMethod("get_standata")
}

#' @export
get_standata.bscmfit <- function(x, ...) {
  stopifnot_(
    !is.null(x$data),
    "Reconstructing Stan data requires the original data.
    Refit the model with {.code save_data = TRUE}."
  )
  Y <- get_stan_y(x)
  Z <- get_stan_Z(x)
  na_y <- which(is.na(Y), arr.ind = TRUE)
  Y[na_y] <- 0
  missing_idx <- cbind(unit = na_y[, 2L], time = na_y[, 1L])
  X_y <- X_z <- NULL
  if (x$setup$has_x) {
    unit_var <- x$setup$unit
    X_mat <- stats::model.matrix(x$formula, data = x$data)
    if (x$setup$has_icpt) {
      X_mat <- X_mat[, -1L, drop = FALSE]
    }
    K <- ncol(X_mat)
    N <- get_N(x)
    J <- get_J(x)
    T_total <- x$setup$T_total
    is_treated <- x$data[[unit_var]] %in% x$setup$treated
    X_y_mat <- X_mat[is_treated, , drop = FALSE]
    X_z_mat <- X_mat[!is_treated, , drop = FALSE]
    X_y <- simplify2array(
      lapply(seq_len(K), \(k) matrix(X_y_mat[, k], N, T_total, TRUE))
    )
    X_z <- simplify2array(
      lapply(seq_len(K), \(k) matrix(X_z_mat[, k], J, T_total, TRUE))
    )
  }

  spline_def <- if (x$setup$has_w) {
    build_spline(
      x$setup$T_total,
      x$setup$T_pre,
      x$setup$spline_df,
      x$setup$spline_type,
      x$setup$noncentered_xi
    )
  }
  standata <- create_standata(
    x$setup,
    x$priors,
    Y,
    Z,
    X_y,
    X_z,
    spline_def,
    x$setup$prior_only
  )
  standata$missing_idx <- missing_idx
  standata$n_missing <- nrow(missing_idx)
  standata
}
