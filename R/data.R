#' Simulated Example Data
#'
#' This simulated data is generated based on a latent factor model with three 
#' factors, 50 time points and 31 units. For the first unit (id 1), the 
#' outcome variable `y` is augmented with a "treatment" `log(1 + 1:10)` for the 
#' last 10 time points.
#' To see exactly how the data was generated, see `data-raw` folder on the 
#' Github repository of the package.
#'
#' @format
#' A data frame with 1,550 rows and 5 columns:
#' \describe{
#'   \item{time}{Time index from ranging from -40 to 9.}
#'   \item{id}{Unit index variable ranging from 1 to 31.}
#'   \item{y}{Outcome variable.}
#'   \item{treatment}{Binary indicator variable where 1 corresponds to treatment.}
#'   \item{tau}{True treatment effect.}
#' }
#'
"simulated_example"