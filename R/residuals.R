#' Posterior residuals of a Bayesian synthetic control model
#'
#' Returns posterior draws or summaries of residuals, defined as the
#' difference between the observed outcome and the posterior expected value.
#'
#' @inheritParams rmse.bscmfit
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param pretreatment_only \[`logical(1)`]\cr If `TRUE` (the default), only
#'   pre-treatment time points are returned, while `FALSE` returns all time
#'   points.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or
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
  probs <- sort_probs(probs)
  test_summary(summary)
  stopifnot_(
    !is.null(object$data),
    "The model fit {.arg object} does not contain the original data. You
    probably used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  d <- fitted(object, summary = FALSE)
  y <- c(get_stan_y(object))
  d <- d |>
    dplyr::mutate(residuals = y - .data$y_mean) |>
    dplyr::select(-"y_mean")
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(.data$residuals, probs)) |>
      dplyr::select(-"residuals", -"variable")
  }
  if (pretreatment_only) {
    d <- d |>
      dplyr::filter(.data[[get_treatment(object)]] == 0)
  }
  d
}
