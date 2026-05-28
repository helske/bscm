#' Extract posterior draws of model parameters as a data frame
#'
#' Returns a `data.frame` representation of the posterior sample of the model
#' parameters.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param parameters \[`character()`]\cr Vector of parameter names.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Ignored.
#' @return A `data.frame` containing model parameters in a wide format.
#' @seealso [as_draws.bscmfit()].
#' @aliases as.data.frame
#' @export
#' @examples
#' head(as.data.frame(fit_single_treated, parameters = c("alpha", "sigma")))
as.data.frame.bscmfit <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  parameters = NULL,
  ...
) {
  all_pars <- setdiff(get_stanfit(x)@model_pars, c("omega_raw", "a"))
  if (is.null(parameters)) {
    parameters <- c("alpha", "beta", "sigma", "omega")
    parameters <- intersect(parameters, all_pars)
  } else {
    stopifnot_(
      checkmate::test_subset(parameters, all_pars),
      "Model does not contain parameters 
      {.val {setdiff(parameters, all_pars)}}."
    )
  }
  as.data.frame(get_stanfit(x), pars = parameters)
}
