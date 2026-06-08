#' Posterior draws of pointwise log-likelihood
#'
#' Returns posterior draws of pointwise log-likelihoods of the treated units per
#' time point. By default (`T_end = NULL`), only pre-treatment observations are
#' included (`T_pre[i]` time points for unit `i`). When `T_end` is supplied as a
#' value of the time variable, all treated units are evaluated up to and
#' including that time point regardless of whether these correspond to pre- or
#' post-treatment times.
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param T_end \[`integer(1)`]\cr A value of the time variable
#' (as it appears in the data) up to which log-likelihoods are computed for all
#' treated units.
#' @param ... Ignored.
#' @return A matrix where rows correspond to posterior draws and columns to
#' time periods. By default, only pre-treatment observations are
#' included for each unit, so output has `sum(T_pre)` columns in total. When
#' `T_end` is supplied, output has `N * T_max` columns where `T_max` is the
#' number of time points from the first period until `T_end`. Output columns are
#' ordered by unit then time: all (pre-treatment) time points for the first
#' treated unit, then the second, and so on.
#' @aliases log_lik
#' @export log_lik
#' @export
#' @examples
#' log_lik(fit_single_treated, T_end = -28) |> head()
log_lik.bscmfit <- function(object, T_end = NULL, ...) {
  T_total <- get_T_total(object)
  T_pre <- get_T_pre(object)
  N <- get_N(object)
  treated <- get_treated(object)
  treatment <- get_treatment(object)
  outcome <- get_outcome(object)
  unit <- get_unit(object)
  mu_draws <- posterior_epred(object)
  sigma_draws <- as.matrix(get_stanfit(object), pars = "sigma")

  stopifnot_(
    !is.null(object$data),
    "Computing log-likelihood requires the original data. Refit the model
      with {.code save_data = TRUE}."
  )

  if (is.null(T_end)) {
    n_ll <- sum(T_pre)
    ll <- matrix(NA_real_, nrow(mu_draws), n_ll)
    idx <- 0L
    for (i in seq_len(N)) {
      id <- treated[i]
      y_pre <- object$data |>
        dplyr::filter(.data[[unit]] == .env$id & .data[[treatment]] == 0) |>
        dplyr::pull(.data[[outcome]])
      mu_offset <- (i - 1L) * T_total
      s <- sigma_draws[, i]
      for (t in seq_along(y_pre)) {
        idx <- idx + 1L
        ll[, idx] <- stats::dnorm(
          y_pre[t],
          mu_draws[, mu_offset + t],
          s,
          log = TRUE
        )
      }
    }
  } else {
    times <- get_times(object)
    T_max <- match(T_end, times)
    stopifnot_(
      !is.na(T_max),
      "Argument {.arg T_end} must be a value of the time variable present in
      the original data."
    )
    n_ll <- T_max * N
    ll <- matrix(NA_real_, nrow(mu_draws), n_ll)
    idx <- 0L
    for (i in seq_len(N)) {
      id <- treated[i]
      y_obs <- object$data |>
        dplyr::filter(.data[[unit]] == .env$id) |>
        dplyr::pull(.data[[outcome]])
      mu_offset <- (i - 1L) * T_total
      s <- sigma_draws[, i]
      for (t in seq_len(T_max)) {
        idx <- idx + 1L
        ll[, idx] <- stats::dnorm(
          y_obs[t],
          mu_draws[, mu_offset + t],
          s,
          log = TRUE
        )
      }
    }
  }
  colnames(ll) <- paste0("log_lik[", seq_len(n_ll), "]")
  ll
}

loglik_at_times <- function(fit, data, t_positions) {
  unit <- get_unit(fit)
  outcome <- get_outcome(fit)
  T_total <- get_T_total(fit)
  N <- get_N(fit)
  treated <- get_treated(fit)
  mu_draws <- posterior_epred(fit)
  sigma_draws <- as.matrix(get_stanfit(fit), pars = "sigma")

  n_t <- length(t_positions)
  ll <- matrix(0, nrow(mu_draws), n_t)
  for (i in seq_len(N)) {
    id <- treated[i]
    y_obs <- data |>
      dplyr::filter(.data[[unit]] == .env$id) |>
      dplyr::pull(.data[[outcome]])
    mu_offset <- (i - 1L) * T_total
    s <- sigma_draws[, i]
    for (j in seq_len(n_t)) {
      ll[, j] <- ll[, j] +
        stats::dnorm(
          y_obs[t_positions[j]],
          mu_draws[, mu_offset + t_positions[j]],
          s,
          log = TRUE
        )
    }
  }
  ll
}
