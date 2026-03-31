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
 
  if (has_predictors(object)) {
    cf <- coef(object, probs = probs)
  } else {
    cf <- NULL
  }
  s <- sigma(object, probs = probs) |> 
    mutate(variable = "Residual SD")
  r2 <- bayes_R2(object, probs = probs) |> 
    mutate(variable = "Bayesian R-squared")
  att <- as_draws(object, parameters = "avg_effect_post") |> 
    summarise_with_probs(probs) |> 
    mutate(variable = "Average treatment effect")
  rmses <- rmse(object, probs = probs, average = TRUE)
  eff <- effective_donors(object, probs = probs, average = TRUE)
  sumr <- bind_rows(cf, s, r2, att, rmses, eff)
  class(sumr) <- c("summary_bscmfit", class(sumr))
  sumr
}
