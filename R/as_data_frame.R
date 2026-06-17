#' Extract posterior draws of model parameters as a data frame
#'
#' Returns a `data.frame` representation of the posterior sample of the model
#' parameters. For samples from posterior predictive distribution, 
#' see [bscm::posterior_epred()] and [bscm::posterior_predict()]. For donor 
#' weights, use [bscm::donor_weights()].
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param parameters \[`character()`]\cr Vector of parameter names. When
#'   `NULL`, (the default), corresponds to a relevant subset of
#'   `c("alpha", "beta", "sigma", "sigma_gamma", "rho")`. Other possible 
#'   choices are `"gamma"` (time-varying regression coefficients) and `"lp__` 
#'   (log-posterior values without constants).
#' @param include \[`logical(1)`]\cr If `TRUE` (the default), output includes
#'   only the variables defined by the argument `parameters`. If `FALSE`, these
#'   variables are excluded from the output. If `NULL`, same as `TRUE` but
#'   variables not present in the model object are silently ignored
#'   (whereas `TRUE` throws an error).
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
  include = TRUE,
  ...
) {
  all_pars <- setdiff(get_stanfit(x)@model_pars, x$setup$excluded_pars)
  if (is.null(parameters)) {
    parameters <- c("alpha", "beta", "sigma", "sigma_gamma", "rho")
    parameters <- intersect(parameters, all_pars)
  } else {
    if (is.null(include)) {
      pars <- sub("\\[.*", "", parameters)
      stopifnot_(
        checkmate::test_subset(unique(pars), all_pars),
        "Model does not contain any of the parameters in {.arg parameters}."
      )
      parameters <- parameters[pars %in% all_pars]
      include <- TRUE
    } else {
      pars <- unique(sub("\\[.*", "", parameters))
      stopifnot_(
        checkmate::test_subset(pars, all_pars),
        "Model does not contain parameters {.val {setdiff(pars, all_pars)}}."
      )
    }
  }
  as.data.frame(get_stanfit(x), pars = parameters)
}
