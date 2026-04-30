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
#' @export posterior_predict
#' @export
posterior_predict.bscmfit <- function(object, ...) {
  as.matrix(get_stanfit(object), pars = "y_rep")
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
#' @export posterior_epred
#' @export
posterior_epred.bscmfit <- function(object, ...) {
  as.matrix(get_stanfit(object), pars = "y_mean")
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
#' @export posterior_linpred
#' @export
posterior_linpred.bscmfit <- function(object, transform = FALSE, ...) {
  as.matrix(get_stanfit(object), pars = "y_mean")
}
