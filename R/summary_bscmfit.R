#' Summarise posterior draws of an estimated Bayesian synthetic control model
#'
#' Generates posterior summary statistics for a Bayesian synthetic control
#' model estimated with [bscm()]. Note that many (more) summaries are available
#' via separate methods, e.g., [donor_weights()].
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A data frame with posterior summaries of
#'   * Estimated intercept terms, time-invariant regression coefficients,
#'     standard deviations of time-varying regression coefficients, and
#'     residual standard deviation parameters \eqn{\sigma}.
#'   * Bayesian \eqn{R^2} values for each treated unit
#'   * Average RMSE for pre- and post-treatment periods
#'   * Average effective number of donors
#'   * Average treatment effect (over post-treatment period and treated units)
#' @aliases summary
#' @export
summary.bscmfit <- function(object, probs = c(0.025, 0.975), ...) {
  probs <- sort_probs(probs)
  cf <- NULL
  if (has_intercept(object) || has_predictors(object)) {
    cf <- coef(object, probs = probs)
    if (length(object$setup$gamma_names) > 0) {
      cf <- cf |>
        dplyr::filter(!startsWith(.data$variable, "gamma")) |>
        dplyr::select(-get_time(object))
    }
  }
  s <- sigma(object, probs = probs) |>
    dplyr::mutate(variable = "Residual SD")
  att <- treatment_effect(
    object,
    type = "average",
    average = FALSE,
    probs = probs
  ) |>
    dplyr::rename(variable = .data$treatment) |>
    dplyr::mutate(variable = c("Pre-treatment effect", "Post-treatment effect"))
  rmses <- rmse(object, average = FALSE, probs = probs) |>
    dplyr::rename(variable = .data$treatment) |>
    dplyr::mutate(variable = c("Pre-treatment RMSE", "Post-treatment RMSE"))
  sumr <- dplyr::bind_rows(cf, s, att, rmses)
  class(sumr) <- c("summary_bscmfit", class(sumr))
  sumr |> dplyr::relocate(.data$variable, .before = 1L)
}
