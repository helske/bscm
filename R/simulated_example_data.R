#' Simulated Example Data
#'
#' This simulated data is generated based on a latent factor model with three 
#' factors, 50 time points and 31 units. For the first unit (`id = "id_1"`), 
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
#'  * id: Unit index variable `"id_1"` to `"id_31"`.
#'  * y: Outcome variable.
#'  * x1: Time-invariant predictor.
#'  * x2: Time-varying predictor.
#'  * treatment: Binary indicator variable where 1 corresponds to the treatment.
#'  * True treatment effect.
#' @name simulated_example
#' @examples
#' head(simulated_example)
NULL