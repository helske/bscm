#' @export
#' @rdname covariate_imbalance
covariate_imbalance <- function(x, ...) {
  UseMethod("covariate_imbalance", x)
}
#' Covariate imbalance of Bayesian synthetic control model
#'
#' For models with covariates, returns and optionally visualizes the
#' covariate imbalances
#' \deqn{\delta_{t} =
#' \sqrt{\frac{1}{K}\sum_{k=1}^K(x_{k,0,t} - \bar x_{k,0,t})^2},}
#' \eqn{t=1,\ldots,T}, where
#' \eqn{\bar x_{k,0,t} = \sum_{j=1}^J \omega_j x_{k, j, t}} and \eqn{x_{k,0,t}}
#' is the value of \eqn{k}th covariate of a treated unit at time t, and
#' similarly for donors \eqn{j=1,\ldots,J}. This is computed separately for
#' each treated unit in case of multiple treated units.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param plot \[`logical(1)`]\cr If `TRUE` (the default), plots the posterior
#' mean and interval of the synthetic covariate distances over time.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#' Default is `c(0.025, 0.975)`. If length of `probs` less than 2, no posterior
#' intervals are drawn, and if length of `probs` is larger than two, the most
#' extreme values are used for the posterior intervals.
#' @param ... Optional arguments passed to [ggplot2::facet_wrap()].
#' @return A `data.frame` of posterior summaries of synthetic covariate
#' distances.
#' @rdname covariate_imbalance
#' @aliases covariate_imbalance
#' @seealso [covariate_adjustment()].
#' @export
#' @examples
#' covariate_imbalance(fit_single_treated, plot = TRUE, probs = c(0.05, 0.95))
#'
covariate_imbalance.bscmfit <- function(
  x,
  plot = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  stopifnot_(
    checkmate::test_flag(plot),
    "Argument {.arg plot} must be a single {.cls logical} value."
  )
  probs <- sort_probs(probs)
  stopifnot_(
    has_predictors(x),
    "The model does not contain any predictors."
  )
  X <- get_Xs(x)
  N <- get_N(x)
  treated <- get_treated(x)
  unit <- get_unit(x)
  delta <- lapply(
    seq_len(N),
    \(i) covariate_imbalance_unit(x, i, X, probs)
  ) |>
    stats::setNames(treated) |>
    dplyr::bind_rows(.id = unit)
  if (plot) {
    time <- get_time(x)
    lwr <- paste0("q", 100 * min(probs))
    upr <- paste0("q", 100 * max(probs))
    ribbon <- if (lwr != upr) {
      geom_ribbon(
        aes(ymin = .data[[lwr]], ymax = .data[[upr]]),
        alpha = 0.25
      )
    }
    facet <- if (N > 1L) {
      facet_wrap(stats::as.formula(paste("~", unit)), ...)
    }
    p <- delta |>
      ggplot(aes(.data[[time]], mean)) +
      ribbon +
      geom_line() +
      labs(x = time, y = "Covariate imbalance") +
      facet +
      theme_bw()
    print(p)
  }
  delta
}

covariate_imbalance_unit <- function(x, unit, X, probs) {
  time <- get_time(x)
  times <- get_times(x)
  T_total <- get_T_total(x)
  J <- get_J(x)
  K <- length(get_predictors(x))

  pars <- paste0("omega[", seq_len(J), ",", unit, "]")
  omega <- posterior::rvar(as_draws(x, pars), with_chains = TRUE)

  X_y <- X$X_y[unit, , , drop = FALSE]
  dim(X_y) <- c(T_total, K)
  delta_tk <- lapply(
    seq_len(K),
    \(k) {
      (X_y[, k] - c(omega %*% X$X_z[,, k]))^2
    }
  )
  sqrt(Reduce(`+`, delta_tk) / K) |>
    posterior::summarise_draws(
      mean,
      ~ posterior::quantile2(.x, probs = probs)
    ) |>
    dplyr::mutate("{time}" := .env$times, .before = "variable") |>
    dplyr::select(-"variable")
}

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
