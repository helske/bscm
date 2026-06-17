#' Return posterior draws from BSCM fit as a draws object 
#'
#' @inheritParams as.data.frame.bscmfit
#' @param inc_warmup \[`logical(1)`]\cr Whether to include warmup draws.
#'   Default is `FALSE`.
#' @return An object of class `draws_array` containing the posterior draws of
#' the specified parameters.
#' @seealso [posterior::as_draws()] for converting the output to other formats,
#' and [posterior::summarise_draws()] for computing posterior summaries of the
#' draws.
#' @aliases as_draws
#' @export as_draws
#' @export
#' @examples
#' head(as_draws(fit_single_treated, parameters = c("alpha", "sigma")))
as_draws.bscmfit <- function(
  x,
  parameters = NULL,
  inc_warmup = FALSE,
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
  posterior::as_draws_array(
    rstan::extract(
      get_stanfit(x),
      pars = parameters,
      permuted = FALSE,
      inc_warmup = inc_warmup,
      include = include
    )
  )
}
