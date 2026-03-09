#' Simulated example data
#'
#' This simulated data is generated based on a latent factor model with three 
#' factors, 50 time points and 31 units. For the first unit (`id = 1`), 
#' a "treatment" \eqn{\tau = log(1 + t)} for the last 10 time points 
#' \eqn{t=0,\ldots,9} is added to an outcome `y`.
#' To see exactly how the data was generated, see `data-raw` folder on the 
#' Github repository of the package.
#'
#' @docType data
#' @keywords datasets
#' @format
#' A data frame with 1,550 rows and 5 columns:
#'  * time: Time index from ranging from -40 to 9.
#'  * id: Unit index variable `ranging from 1 to 31.
#'  * y: Outcome variable.
#'  * x: Time-varying predictor.
#'  * treatment: Binary indicator variable where 1 corresponds to the treatment.
#'  * tau: True treatment effect.
#'  * beta: True time-varying regression coefficient.
#'  * psi1: First latent factor.
#'  * psi2: Second latent factor.
#'  * lambda1: Loadings of the first factor.
#'  * lambda2: Loadings of the second factor.
#' @name simulated_data
#' @examples
#' head(simulated_data)
NULL