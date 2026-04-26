#' Extract residual standard deviations of a Bayesian synthetic control
#' model
#'
#' @inheritParams rmse.bscmfit
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the residual standard
#'   deviations.
#' @aliases sigma
#' @export
sigma.bscmfit <- function(object, probs = c(0.025, 0.975), ...) {
  
  probs <- sort_probs(probs)
  
  d <- as_draws(object, "sigma")
  N <- get_N(object)
  treated <- get_treated(object)
  unit <- get_unit(object)
  as_draws(object, "sigma") |> 
    summarise_with_probs(probs) |>
    mutate("{unit}" := treated, .before = 1L) |> 
    mutate(variable = "sigma")
}
