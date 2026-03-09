#' Extract posterior draws of model parameters as a data frame
#'
#' Returns a `data.frame` representation of the posterior sample of the model
#' parameters.
#' 
#' @export
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param parameters \[`character()`]\cr Vector of parameter names.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Ignored.
#' @return A `data.frame` containing model parameters in a wide format.
#' @seealso [as_draws.bscmfit()].
as.data.frame.bscmfit <- function(x, row.names = NULL, optional = FALSE, 
                                  parameters = NULL, ...) {
  all_pars <- setdiff(x$stanfit@model_pars, c("omega_raw", "a"))
  if (is.null(parameters)) {
    parameters <- c("alpha", "beta", "sigma", "omega", "tau")
    parameters <- intersect(parameters, all_pars)
  } else {
    stopifnot_(
      checkmate::test_subset(parameters, all_pars),
      "Model does not contain parameters 
      {.val {setdiff(parameters, all_pars)}}."
    )
  }
  as.data.frame(x$stanfit, pars = parameters)
}