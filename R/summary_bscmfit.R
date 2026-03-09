#' Summarize posterior draws of an estimated Bayesian synthetic control model
#'
#' Generates posterior summary statistics for a Bayesian synthetic control 
#' model estimated with [bscm()]. The output is organized as a list containing 
#' data frames for different types of summaries.
#'
#' @details
#' By default, the output contains various posterior summary statistics. You can
#' select a subset of these via the `include` argument with the following 
#' accepted values:
#' 
#' * `"effects"`: Treatment effects over time (observed series - synthetic control).
#' * `"synthetic"`: Synthetic control values over time.
#' * `"cumulative_effects"`: Cumulative average treatment effects in the 
#'   post-treatment period.
#' * `"average_effects"`: Average treatment effects in pre- and post-treatment 
#'   periods.
#' * `"relative_change"`: Relative change from last pre-treatment synthetic 
#'    control value in the post-treatment period.
#' * `"parameters"`: Model parameters excluding the weights, e.g., 
#'   intercept (`alpha`), standard deviation of noise term (`sigma`); and in a 
#'   case of covariates, regression coefficients (`beta`).
#' * `"weights"`: Donor weights used in the synthetic control.
#' * `"effective_donors"`: Effective number of donor units based on weight 
#'   concentration \eqn{1/\sum(w^2)}.
#' * `"RMSE"`: Root mean squared errors for pre- and post-treatment periods, and 
#'   their ratio.
#' * `"R2"`: Bayesian \eqn{R^2} measure of model fit.
#'
#' Scalar quantities (`average_effects`, `RMSE`, `R2`, `effective_donors`) are
#' combined into a single `misc` data frame when multiple are requested.
#' Non-scalar quantities are returned as separate data frames indexed by time 
#' or donor unit.
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param include \[`character()` or `NULL`]\cr Character vector defining which 
#'   posterior summaries to include in the output. See details for accepted 
#'   values. If `NULL` (the default), all available summaries are included.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A named list of data frames containing the requested posterior 
#'   summaries. The `misc` element combines all scalar summaries with a 
#'   `variable` column indicating the quantity.
#' @export
summary.bscmfit <- function(object, include = NULL, 
                            probs = c(0.025, 0.975), ...) {
  
  stopifnot_(
    checkmate::test_numeric(
      probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values between
     0 and 1."
  )
  all_output <- c(
    "effects", "synthetic", "weights", "cumulative_effects", "average_effects",
    "relative_change", "parameters", "RMSE", "R2", "effective_donors"
  )
  if (is.null(include)) {
    include <- all_output
  } else {
    matched <- pmatch(include, all_output)
    stopifnot_(
      !anyNA(matched),
      c(
        "Argument {.arg include} contains invalid values.",
        "i" = "Accepted values are: {.val {all_output}}."
      )
    )
    include <- all_output[matched]
  }
  
  # Extract setup info
  setup <- object$setup
  times <- setup$times
  donors <- setup$donors
  T_pre <- setup$T_pre
  T_total <- setup$T_total
  time <- setup$time
  unit <- setup$unit
  has_predictors <- length(setup$predictors) > 0L
  
  # Get draws from Stan fit
  all_pars <- object$stanfit@model_pars
  pars <- intersect(
    c("alpha", "beta", "sigma", "tau"), all_pars
  )
  vars <- c(
    if ("effects" %in% include) "effect",
    if ("synthetic" %in% include) "synthetic_y",
    if ("cumulative_effects" %in% include) "avg_effect_post_cumulative",
    if ("weights" %in% include) "omega",
    if ("relative_change" %in% include) "relative_change",
    if ("average_effects" %in% include) c("avg_effect_pre", "avg_effect_post"),
    if ("parameters" %in% include) pars,
    if ("RMSE" %in% include) c("RMSE_pre", "RMSE_post", "RMSE_ratio"),
    if ("R2" %in% include) "R2",
    if ("effective_donors" %in% include) "effective_donors"
  )
  draws <- as_draws(object, vars)
  
  out <- list()
  
  if ("effects" %in% include) {
    out$effects <- draws |> 
      subset_draws("effect") |> 
      summarize_with_probs(probs) |> 
      mutate("{time}" := .env$times, .before = 1L) |> 
      select(-"variable")
  }
  if ("synthetic" %in% include) {
    out$synthetic <- draws |> 
      subset_draws("synthetic_y") |> 
      summarize_with_probs(probs) |> 
      mutate("{time}" := .env$times, .before = 1L) |> 
      select(-"variable")
  }
  if ("cumulative_effects" %in% include) {
    post_times <- times[(T_pre + 1L):T_total]
    out$cumulative_effects <- draws |> 
      subset_draws("avg_effect_post_cumulative") |> 
      summarize_with_probs(probs) |> 
      mutate("{time}" := .env$post_times, .before = 1L) |> 
      select(-"variable")
  }
  if ("weights" %in% include) {
    out$weights <- draws |> 
      subset_draws("omega") |> 
      summarize_with_probs(probs) |> 
      mutate("{unit}" := .env$donors, .before = 1L) |> 
      select(-"variable")
  }
  if ("relative_change" %in% include) {
    post_times <- times[(T_pre + 1L):T_total]
    out$relative_change <- draws |> 
      subset_draws("relative_change") |> 
      summarize_with_probs(probs) |> 
      mutate("{time}" := .env$post_times, .before = 1L) |> 
      select(-"variable")
  }
  if ("parameters" %in% include) {
    out$parameters <- draws |> 
      subset_draws(pars) |> 
      summarize_with_probs(probs)
  }
  vars <- setdiff(
    vars, 
    c("effect", "synthetic_y", "avg_effect_post_cumulative", "omega", 
      "relative_change", pars)
  )
  if (length(vars) > 0) {
    out$misc <- draws |> 
      subset_draws(vars) |> 
      summarize_with_probs(probs)
  }
  class(out) <- "bscmfit_summary"
  out
}