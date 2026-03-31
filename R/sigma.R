#' Extract residual standard deviations of a Bayesian synthetic control
#' model
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the residual standard
#'   deviations.
#' @aliases sigma
#' @export
sigma.bscmfit <- function(object, probs = c(0.025, 0.975), ...) {
  
  test_probs(probs)
  
  d <- as_draws(object, "sigma")
  N <- get_N(object)
  treated <- get_treated(object)
  unit <- get_unit(object)
  as_draws(object, "sigma") |> 
    summarise_with_probs(probs) |>
    mutate("{unit}" := treated, .before = 1L) |> 
    mutate(variable = "sigma")
}
