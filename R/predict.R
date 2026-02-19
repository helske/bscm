#' Posterior Predictive Distribution of Bayesian Synthetic Control Model
#' 
#' Returns draws from the posterior predictive distribution of the Bayesian 
#' synthetic control. Note that this function does not simulate new 
#' realizations, but rather returns the posterior draws computed during the 
#' model estimation.
#' 
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A matrix of posterior predictive draws of the synthetic control, 
#' where rows correspond to posterior draws and columns to time periods.
#' @export
predict.bscmfit <- function(object, ...) {
  as.matrix(object$stanfit, pars = "synthetic_y")
}
