#' Simulated example data with single treated unit
#'
#' This simulated data is generated based on a latent factor model with two
#' factors, 40 time points and 31 units. For the first unit (`id = 1`),
#' a "treatment" \eqn{\tau = 1 + t} for the last 10 time points
#' \eqn{t=0,\ldots,9} is added to an outcome `y`.
#' To see exactly how the data was generated, see `data-raw` folder on the
#' Github repository of the package.
#'
#' @docType data
#' @keywords datasets
#' @format
#' A data frame with 1,550 rows and 5 columns:
#'  * time: Time index from ranging from -30 to 9.
#'  * id: Unit index variable `ranging from 1 to 31.
#'  * y: Outcome variable.
#'  * x: Time-varying predictor.
#'  * treatment: Binary indicator variable where 1 corresponds to the treatment.
#'  * alpha: Unit-specific intercepts.
#'  * psi1: First latent factor.
#'  * psi2: Second latent factor.
#'  * lambda1: Loadings of the first factor.
#'  * lambda2: Loadings of the second factor.
#' @name single_treated
#' @examples
#' head(single_treated)
#' single_treated |>
#'   ggplot2::ggplot(ggplot2::aes(time, y, group = id)) +
#'   ggplot2::geom_line(alpha = 0.3) +
#'   ggplot2::geom_line(
#'     data = single_treated |> dplyr::filter(id == 1), colour = "tomato"
#'   ) +
#'   ggplot2::theme_bw()
NULL
