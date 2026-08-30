#' @export
#' @rdname rmse
rmse <- function(x, ...) {
  UseMethod("rmse", x)
}
#' Extract root mean squared errors of a Bayesian synthetic control model
#'
#' Returns posterior draws (or summaries) of root mean squared errors (RMSEs)
#' for the pre-treatment and post-treatment periods.
#'
#' @inheritParams bscm_postprocessing
#' @param average \[`logical(1)`]\cr If `TRUE`, returns the
#' average RMSEs over treated units in case of multiple treated units.
#' If `FALSE` (the default), unit-specific values are returned.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`).
#' @rdname rmse
#' @aliases rmse
#' @export
#' @examples
#' rmse(fit_single_treated, probs = c(0.01, 0.5, 0.8))
rmse.bscmfit <- function(
  x,
  average = FALSE,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  check_flag(average, "average")
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  check_has_data(x, "x")
  d <- rmse_draws(x, average)
  if (summary) {
    d <- summarise_column(d, "rmse", probs)
  }
  d
}
#' @export
#' @rdname rmse_ratio
rmse_ratio <- function(x, ...) {
  UseMethod("rmse_ratio", x)
}
#' Extract the ratio of post- to pre-treatment RMSE
#'
#' Returns posterior draws (or summaries) of the ratio of the
#' post-treatment RMSE to the pre-treatment RMSE.
#'
#' @inheritParams bscm_postprocessing
#' @param average \[`logical(1)`]\cr If `TRUE`, returns the
#' average ratio over treated units in case of multiple treated units.
#' If `FALSE` (the default), unit-specific values are returned.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @export
#' @rdname rmse_ratio
#' @examples
#' rmse_ratio(fit_single_treated, probs = c(0.01, 0.5, 0.8))
rmse_ratio.bscmfit <- function(
  x,
  average = FALSE,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  check_flag(average, "average")
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  unit <- get_unit(x)
  treatment <- get_treatment(x)
  d <- rmse_draws(x) |>
    dplyr::summarise(
      ratio = .data$rmse[.data[[treatment]] == 1] /
        .data$rmse[.data[[treatment]] == 0],
      .by = dplyr::all_of(unit)
    )
  if (average && get_N(x) > 1L) {
    d <- d |>
      dplyr::summarise(ratio = posterior::rvar_mean(.data$ratio))
  }
  if (summary) {
    d <- summarise_column(d, "ratio", probs)
  }
  d
}

#' Posterior draws of the pre- and post-treatment RMSEs of a BSCM
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param average \[`logical(1)`]\cr If `TRUE`, the RMSEs are averaged over the
#'   treated units within each posterior draw.
#' @noRd
rmse_draws <- function(x, average = FALSE) {
  treatment <- get_treatment(x)
  d <- effects_draws(x) |>
    dplyr::summarise(
      rmse = sqrt(posterior::rvar_mean(.data$effect^2)),
      .by = dplyr::all_of(c(get_unit(x), treatment))
    )
  if (average && get_N(x) > 1L) {
    d <- d |>
      dplyr::summarise(
        rmse = posterior::rvar_mean(.data$rmse),
        .by = dplyr::all_of(treatment)
      )
  }
  d
}
