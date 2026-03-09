#' Return posterior draws from BSCM fit
#' 
#' @aliases as_draws
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param parameters \[`character()`]\cr Vector of parameter names.
#' @param inc_warmup \[`logical(1)`]\cr Whether to include warmup draws. 
#' Default is `FALSE`.
#' @param include \[`logical(1)`]\cr If `TRUE` (the default), output includes 
#' only the variables defined by the argument `parameters`. If `FALSE`, these 
#' variables are excluded from the output.
#' @param ... Ignored.
#' @return An object of class `draws_array` containing the posterior draws of 
#' the specified parameters.
#' @seealso [posterior::as_draws()] for converting the output to other formats.
#' @export
#' @export as_draws
as_draws.bscmfit <- function(x, parameters, inc_warmup = FALSE, include = TRUE, 
                             ...) {
  posterior::as_draws_array(
    rstan::extract(
      x$stanfit,
      pars = parameters,
      permuted = FALSE,
      inc_warmup = inc_warmup,
      include = include
    )
  )
}
