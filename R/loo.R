#' Approximate leave-One-Out (LOO) cross-validation for Bayesian synthetic control model
#' 
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param r_eff \[`logical(1)`]\cr If `TRUE` (the default), [loo::loo()] 
#' computes more accurate Monte Carlo error estimates at the cost of 
#' increased computation time.
#' @param ... Additional arguments to [rstan::loo()].
#' @return An output from [rstan::loo()].
#' @references Vehtari A, Gelman A, and Gabry J (2017).
#' Practical Bayesian model evaluation using leave-one-out cross-validation and
#' WAIC. *Statistics and Computing*. 27(5), 1413--1432.
#' @aliases loo
#' @export loo
#' @export
loo.bscmfit <- function(x, r_eff = TRUE,...) {
  rstan::loo(get_stanfit(x), r_eff = r_eff, ...)
}
