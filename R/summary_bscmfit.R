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
#'   * Estimated intercept terms, regression coefficients and residual SDs
#'   * Bayesian \eqn{R^2} values for each treated unit
#'   * Average RMSE for pre- and post-treatment periods and their ratio
#'   * Average effective number of donors
#'   * Average treatment effect (over post-treatment period and treated units)
#' @aliases summary
#' @export
summary.bscmfit <- function(object, probs = c(0.025, 0.975), ...) {
  
  test_probs(probs)
  cf <- sigma_gamma <- NULL
  if (has_predictors(object)) {
    cf <- coef(object, probs = probs, type = "beta")$beta |> 
      mutate("{get_unit(object)}" := NA, .before = 1L)
    if (has_tv_coefs(object)) {
      sigma_gamma <- as_draws(object, "sigma_gamma") |> 
        summarise_with_probs(probs) |> 
        mutate(variable = paste0("SD sigma_gamma_", object$setup$gamma_names))
    }
  }
  s <- sigma(object, probs = probs) |> 
    mutate(variable = "Residual SD (sigma)")
  r2 <- bayes_R2(object, probs = probs) |> 
    mutate(variable = "Bayesian R-squared")
  att <- treatment_effect(
    object, type = "average", average = TRUE, probs = probs
  )
  rmses <- rmse(object, average = TRUE, probs = probs)
  eff <- effective_donors(object, average = TRUE, probs = probs)
  sumr <- bind_rows(cf, sigma_gamma, s, r2, att, rmses, eff)
  class(sumr) <- c("summary_bscmfit", class(sumr))
  sumr
}
