#' Extract the `stanfit` object from the `bscmfit` object
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
    is.list(x),
    "Input {.arg x} is not a list-like object."
  )
  x$stanfit
}
#' @export
#' @rdname get_stanfit
get_stanfit.bscmfit <- function(x, ...) {
  x$stanfit
}

#' Return the number of posterior draws of a `bscmfit` object
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
#' Return the number of chains of `bscmfit` object
#' @aliases nchains
#' @export
#' @export nchains
#' @param x \[`bscmfit`]\cr The model fit object.
#' @return Number of Markov chains used in sampling as a single `integer` value.
nchains.bscmfit <- function(x) {
  as.integer(x$stanfit@sim$chains)
}
