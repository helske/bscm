#' Extract residual standard deviations of a Bayesian synthetic control
#' model
#'
#' @inheritParams bscm_postprocessing
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
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
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  unit <- get_unit(object)
  d <- dplyr::tibble(
    "{unit}" := get_treated(object),
    sigma = rvars_of(object, "sigma")
  )
  if (summary) {
    d <- summarise_column(d, "sigma", probs)
  }
  d
}
