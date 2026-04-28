#' Simulated example data with multiple treated units
#'
#' This simulated data is generated based on a latent factor model with two
#' factors, 40 time points and 53 units. For the first three units,
#' a "treatment" \eqn{\tau = 1 + t} for the last 5 time points
#' \eqn{t=0,\ldots,4} is added to an outcome `y`.
#' To see exactly how the data was generated, see `data-raw` folder on the
#' Github repository of the package.
#'
#' @docType data
#' @keywords datasets
#' @format
#' A data frame with 1,550 rows and 5 columns:
#'  * time: Time index from ranging from -30 to 9.
#'  * id: Unit index variable `ranging from 1 to 33.
#'  * y: Outcome variable.
#'  * x: Time-varying predictor.
#'  * z: Time-varying predictor.
#'  * treatment: Binary indicator variable where 1 corresponds to the treatment.
#'  * tau: True treatment effect.
#'  * psi1: First latent factor.
#'  * psi2: Second latent factor.
#'  * lambda1: Loadings of the first factor.
#'  * lambda2: Loadings of the second factor.
#' @name multiple_treated
#' @examples
#' head(multiple_treated)
NULL
