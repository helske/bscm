#' Expected values of posterior predictive distribution of Bayesian synthetic
#' control model
#'
#' Returns posterior draws (or summaries) of expected values `y_mean` of the
#' posterior predictive distribution of the model.
#'
#' @inheritParams bscm_postprocessing
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @aliases fitted
#' @export
#' @examples
#' fitted(fit_single_treated) |> head()
fitted.bscmfit <- function(
  object,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  check_has_data(object, "object")
  d <- fitted_draws(object)
  if (summary) {
    d <- summarise_column(d, "y_mean", probs)
  }
  d
}

#' Posterior draws of the expected values of a BSCM
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @noRd
fitted_draws <- function(x) {
  y_mean <- gq_to_rvar(x, x$y_mean)
  treated_grid(x) |>
    dplyr::mutate(y_mean = c(.env$y_mean))
}
