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
#'   * If part of the estimated model, posterior summaries of intercept terms,
#'     time-invariant regression coefficients, standard deviations of
#'     time-varying regression coefficients, residual standard deviations,
#'     and autoregressive coefficients of residuals.
#'   * Average RMSE for pre- and post-treatment periods
#'   * Average treatment effect (over post-treatment period and treated units)
#' @aliases summary
#' @export
summary.bscmfit <- function(object, probs = c(0.025, 0.975), ...) {
  probs <- sort_probs(probs)
  cf <- rho <- NULL
  ar1_error <- object$setup$error == "ar1"
  if (has_intercept(object) || has_predictors(object) || ar1_error) {
    cf <- coef(
      object,
      probs = probs,
      type = c("alpha", "beta", "sigma_gamma", "rho")
    )
  }
  s <- sigma(object, probs = probs) |>
    dplyr::mutate(variable = "Residual SD")
  treatment <- get_treatment(object)
  att <- treatment_effect(
    object,
    type = "average",
    average = FALSE,
    probs = probs
  ) |>
    dplyr::rename(variable = treatment) |>
    dplyr::mutate(
      variable = dplyr::if_else(
        .data$variable == 0,
        "Pre-treatment effect",
        "Post-treatment effect"
      )
    )
  rmses <- rmse(object, average = FALSE, probs = probs) |>
    dplyr::rename(variable = treatment) |>
    dplyr::mutate(
      variable = dplyr::if_else(
        .data$variable == 0,
        "Pre-treatment RMSE",
        "Post-treatment RMSE"
      )
    )
  sumr <- dplyr::bind_rows(cf, s, att, rmses)
  class(sumr) <- c("summary_bscmfit", class(sumr))
  sumr |> dplyr::relocate("variable", .before = 1L)
}
