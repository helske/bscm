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
#' @return A data frame with various posterior summaries of the model.
#' @aliases summary
#' @export
summary.bscmfit <- function(object, probs = c(0.025, 0.975), ...) {
  probs <- sort_probs(probs)
  cf <- rho <- NULL
  coefs <- c(
    if (has_intercept(object)) "alpha",
    if (has_predictors(object)) "beta",
    if (has_tv_coefs(object)) "kappa",
    if (has_ar1_error(object)) "rho"
  )
  if (length(coefs) > 0) {
    cf <- coef(
      object,
      probs = probs,
      type = coefs
    )
  }
  s <- sigma(object, probs = probs) |>
    dplyr::mutate(variable = "Residual SD")
  r2 <- bayes_R2(object, probs = probs) |>
    dplyr::mutate(variable = "Bayesian R2")
  eff <- effective_donors(object, probs = probs) |>
    dplyr::mutate(variable = "Effective number of donors")
  treatment <- get_treatment(object)
  att <- average_treatment_effect(
    object,
    average = FALSE,
    probs = probs
  ) |>
    dplyr::rename(variable = treatment) |>
    dplyr::mutate(
      variable = dplyr::if_else(
        .data$variable == 0,
        "Average pre-treatment effect",
        "Average post-treatment effect"
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
  sumr <- dplyr::bind_rows(cf, s, r2, eff, att, rmses) |>
    dplyr::relocate("variable", .before = 1L)
  class(sumr) <- c("summary_bscmfit", class(sumr))
  attr(sumr, "N") <- get_N(object)
  attr(sumr, "unit") <- get_unit(object)
  sumr
}
