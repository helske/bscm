#' Extract the `stanfit` Object from the `bscmfit` Object
#' 
#' This function returns the output (`stanfit` object) from [rstan::sampling()].
#' 
#' @export
#' @rdname get_stanfit
#' @param x \[`bscmfit`]\cr The output returned by the [bscm()].
#' @param ... Ignored.
#' @return Object of class`stanfit`.
get_stanfit <- function(x, ...) {
  UseMethod("get_stanfit", x)
}
#' @export
#' @rdname get_stanfit
get_stanfit.default <- function(x, ...) {
  stopifnot_(
    !is.list(x),
    "Input {.arg x} is not a list-like object."
  )
  x$stanfit
}
#' @export
#' @rdname get_stanfit
get_stanfit.bscmfit <- function(x, ...) {
  x$stanfit
}

#' Return the Number of Posterior Draws of a `bscmfit` Object
#'
#' @aliases ndraws
#' @export
#' @export ndraws
#' @param x \[`bscmfit`]\cr The model fit object.
#' @return Number of posterior draws as a single `integer` value.
ndraws.bscmfit <- function(x) {
  as.integer((x$stanfit@sim$n_save[1L] - x$stanfit@sim$warmup2[1L]) * 
               x$stanfit@sim$chains)
}
#' Return the Number of Chains of `bscmfit` Object
#' @aliases nchains
#' @export
#' @export nchains
#' @param x \[`bscmfit`]\cr The model fit object.
#' @return Number of Markov chains used in sampling as a single `integer` value.
nchains.bscmfit <- function(x) {
  as.integer(x$stanfit@sim$chains)
}
#' Return Posterior Draws from BSCM fit
#' 
#' @aliases as_draws
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param variable \[`character()`]\cr Vector of parameter names.
#' @param inc_warmup \[`logical(1)`]\cr Whether to include warmup draws. 
#' Default is `FALSE`.
#' @param include \[`logical(1)`]\cr If `TRUE` (the default), output includes 
#' only the variables defined by the argument `variable`. If `FALSE`, these 
#' variables are excluded from the output.
#' @param ... Ignored.
#' @return An object of class `draws_array` containing the posterior draws of 
#' the specified parameters.
#' @seealso [posterior::as_draws()] for converting the output to other formats.
#' @export
#' @export as_draws
as_draws.bscmfit <- function(x, variable, inc_warmup = FALSE, include = TRUE, 
                             ...) {
  posterior::as_draws_array(
    rstan::extract(
      x$stanfit,
      pars = variable,
      permuted = FALSE,
      inc_warmup = inc_warmup,
      include = include
    )
  )
}
