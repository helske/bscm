#' @export
#' @rdname variance_components
variance_components <- function(x, ...) {
  UseMethod("variance_components", x)
}
#' Decomposition of predictive variance
#'
#' Decomposes the posterior predictive variance of the outcome into
#' uncertainty in the expected outcome and observation noise, separately
#' for each treated unit and time point.
#' 
#' @inheritParams as.data.frame.bscmfit
#' @return A `tibble` with the unit, time, and treatment indicators, and
#'   \describe{
#'     \item{`var_epred`}{
#'       Posterior variance of the expected outcome. This measures
#'       uncertainty in the expected outcome across posterior draws.
#'     }
#'     \item{`var_error`}{
#'       Posterior expected observation variance, \eqn{E(\sigma^2 \mid y)}.
#'     }
#'     \item{`var_total`}{
#'       Posterior predictive variance, equal to `var_epred + var_error`.
#'     }
#'     \item{`var_prop`}{
#'       Proportion of posterior predictive variance attributable to
#'       uncertainty in the expected outcome, equal to
#'       `var_epred / var_total`.
#'     }
#'   }
#' @examples
#' \dontrun{
#' fit <- bscm(y ~ treatment, unit = id, time = time, data = d)
#' variance_components(fit_single_treated)
#' }
#' @rdname variance_components
#' @export
variance_components.bscmfit <- function(x, ...) {
  y_mean <- gq_to_rvar(x, x$y_mean)
  sigma <- rvars_of(x, "sigma")
  T_total <- get_T_total(x)

  treated_grid(x) |>
    dplyr::mutate(
      var_epred = posterior::var(c(.env$y_mean)),
      var_error = rep(posterior::E(.env$sigma^2), each = .env$T_total),
      var_total = .data$var_epred + .data$var_error,
      var_prop = .data$var_epred / .data$var_total
    )
}
