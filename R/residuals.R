#' Posterior residuals of a Bayesian synthetic control model
#'
#' Returns posterior draws or summaries of residuals, defined as the
#' difference between the observed outcome and the posterior expected value.
#'
#' @inheritParams bscm_postprocessing
#' @param pretreatment_only \[`logical(1)`]\cr If `TRUE` (the default), only
#'   pre-treatment time points are returned, while `FALSE` returns all time
#'   points.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @aliases residuals
#' @export
#' @examples
#' residuals(fit_single_treated) |> head()
residuals.bscmfit <- function(
  object,
  summary = TRUE,
  probs = c(0.025, 0.975),
  pretreatment_only = TRUE,
  ...
) {
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  check_has_data(object, "object")
  d <- residuals_draws(object, pretreatment_only)
  if (summary) {
    d <- summarise_column(d, "residuals", probs)
  }
  d
}

#' Posterior draws of the residuals of a BSCM
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param pretreatment_only \[`logical(1)`]\cr If `TRUE` (the default), only
#'   pre-treatment time points are returned.
#' @noRd
residuals_draws <- function(x, pretreatment_only = TRUE) {
  y <- c(get_stan_y(x))
  d <- fitted_draws(x) |>
    dplyr::mutate(residuals = .env$y - .data$y_mean) |>
    dplyr::select(-"y_mean")
  if (pretreatment_only) {
    d <- d |>
      dplyr::filter(.data[[get_treatment(x)]] == 0)
  }
  d
}
