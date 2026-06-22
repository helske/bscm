#' Simulated example data with multiple treated units
#'
#' This simulated data is generated based on a latent factor model with two
#' factors, 40 time points and 53 units. For the first three units,
#' a "treatment" \eqn{\tau = 2} for the last 5 time points
#' \eqn{t=0,\ldots,4} is added to an outcome `y`. There are two time-varying
#' covariates `x` and `z`, with latter one having time-varying effect on `y`.
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
#' @name multiple_treated
#' @examples
#' head(multiple_treated)
#' multiple_treated |>
#'   dplyr::mutate(treated = any(treatment), .by = id) |>
#'   ggplot2::ggplot(ggplot2::aes(time, y)) +
#'   ggplot2::geom_line(
#'     ggplot2::aes(group = id, colour = treated, alpha = treated)
#'   ) +
#'   ggplot2::scale_alpha_manual(values = c(0.3, 1)) +
#'   ggplot2::theme_bw()
NULL
