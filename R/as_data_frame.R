#' Extract Posterior Draws of Model Parameters as a Data Frame
#'
#' Returns a `data.frame` representation of the posterior sample of the model
#' parameters.
#' 
#' @export
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Ignored.
#' @return A `data.frame` containing model parameters in a wide format.
#' @seealso [as_draws.bscmfit()].
as.data.frame.bscmfit <- function(x, row.names = NULL, optional = FALSE, ...) {
  pars <- c(
    "alpha", "beta", "sigma", "omega", "alpha_z", "sigma_z", "sigma_delta",
    "mu", "gamma"
  )
  pars <- intersect(pars, x$stanfit@model_pars)
  as.data.frame(x$stanfit, pars = pars)
}