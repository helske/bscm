#' Posterior Predictive Distribution of Bayesian Synthetic Control Model
#' 
#' Returns draws from the posterior predictive distribution of the Bayesian 
#' synthetic control. Note that this function does not simulate new 
#' realizations from this distribution, but rather returns the posterior draws 
#' computed during the model estimation.
#' 
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A matrix of posterior predictive draws of the synthetic control, 
#' where rows correspond to posterior draws and columns to time periods.
#' @aliases posterior_predict
#' @export
#' @export posterior_predict
posterior_predict.bscmfit <- function(object, ...) {
  as.matrix(object$stanfit, pars = "synthetic_y")
}
#' Posterior Draws of the Expected Predictive Distribution
#'
#' Returns draws from the posterior distribution of the expected value of
#' the Bayesian synthetic control model.
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A matrix where rows correspond to posterior draws and columns to
#'   time periods.
#' @aliases posterior_epred
#' @export
#' @export posterior_epred
posterior_epred.bscmfit <- function(object, ...) {
  as.matrix(object$stanfit, pars = "synthetic_mean")
}
#' Posterior Draws of the Linear Predictor
#'
#' Returns draws from the posterior distribution of the linear predictor of
#' the Bayesian synthetic control model. Since the model uses a Gaussian
#' likelihood with identity link, this is equivalent to
#' [posterior_epred.bscmfit()].
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param transform \[`logical(1)`]\cr Ignored.
#' @param ... Ignored.
#' @return A matrix where rows correspond to posterior draws and columns to
#' time periods.
#' @aliases posterior_linpred
#' @export
#' @export posterior_linpred
posterior_linpred.bscmfit <- function(object, transform = FALSE, ...) {
  as.matrix(object$stanfit, pars = "synthetic_mean")
}
#' Posterior draws of pointwise log-likelihood
#'
#' Returns draws of log-likelihood values of the treated unit in the 
#' time points corresponding to the pre-treatment period.
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A matrix where rows correspond to posterior draws and columns to
#' time periods.
#' @aliases log_lik
#' @export
#' @export log_lik
log_lik <- function(object, ...) {
  as.matrix(object$stanfit, pars = "log_lik")
}
