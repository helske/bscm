#' Extract residual standard deviations of a Bayesian synthetic control
#' model
#'
#' @inheritParams rmse.bscmfit
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @aliases sigma
#' @export
#' @examples
#' sigma(fit_single_treated)
sigma.bscmfit <- function(
  object,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  test_summary(summary)
  probs <- sort_probs(probs)
  treated <- get_treated(object)
  unit <- get_unit(object)
  d <- tibble(
    "{unit}" := treated,
    sigma = posterior::as_draws_rvars(as_draws(object, "sigma"))$sigma
  )
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(.data$sigma, probs)) |>
      dplyr::select(-"sigma", -"variable")
  }
  if (length(treated) == 1) {
    d <- d |> dplyr::select(-dplyr::all_of(unit))
  }
  d
}
