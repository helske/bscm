#' Simulated example data with single treated unit
#'
#' This simulated data is generated based on a latent factor model with three
#' factors, 50 time points and 51 units. For the first unit (`id = 1`),
#' a "treatment" \eqn{\tau = 1 + t} for the last 10 time points
#' \eqn{t=0,\ldots,9} is added to an outcome `y`.
#' To see exactly how the data was generated, see `data-raw` folder on the
#' Github repository of the package.
#' @docType data
#' @keywords datasets
#' @format
#' A data frame with 2550 rows and 5 columns:
#'  * time: Time index from ranging from -40 to 9.
#'  * id: Unit index variable `ranging from 1 to 31.
#'  * y: Outcome variable.
#'  * x: Time-varying predictor.
#'  * treatment: Binary indicator variable where 1 corresponds to the treatment.
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
#' Example bscmfit object
#'
#' The object `fit_single_treated` contains a Bayesian synthetic
#' control model estimated as
#' \preformatted{
#' fit <- bscm(
#'   formula = y ~ x, data = single_treated, treatment = "treatment",
#'   chains = 2, cores = 1, refresh = 0, iter = 2000, warmup = 1000,
#'   control = list(adapt_delta = 0.8)
#' )
#' }
#' @docType data
#' @format A `bscmfit` object.
#' @name fit_single_treated
#' @examples
#' fit_single_treated
#' plot(fit_single_treated)
NULL
