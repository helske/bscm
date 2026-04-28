#' Approximate leave-One-Out (LOO) cross-validation for Bayesian synthetic control model
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param r_eff \[`logical(1)`]\cr If `TRUE` (the default), [loo::loo()]
#' computes more accurate Monte Carlo error estimates at the cost of
#' increased computation time.
#' @param ... Additional arguments to [loo::loo()].
#' @return An output from [loo::loo()].
#' @references Vehtari A, Gelman A, and Gabry J (2017).
#' Practical Bayesian model evaluation using leave-one-out cross-validation and
#' WAIC. *Statistics and Computing*. 27(5), 1413--1432.
#' @aliases loo
#' @export loo
#' @export
loo.bscmfit <- function(x, r_eff = TRUE, ...) {
  ll <- log_lik(x)
  nc <- nchains(x)
  cid <- rep(seq_len(nc), each = ndraws(x) / nc)
  r_eff_val <- if (r_eff) loo::relative_eff(exp(ll), chain_id = cid) else 1
  loo::loo(ll, r_eff = r_eff_val, ...)
}
