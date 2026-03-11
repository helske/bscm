#' Summarize posterior draws of an estimated Bayesian synthetic control model
#'
#' Generates posterior summary statistics for a Bayesian synthetic control 
#' model estimated with [bscm()]. Note that most of the summaries are also 
#' available via separate methods, e.g., [donor_weights()].
#'
#' @details
#' By default, the output contains various posterior summary statistics. You can
#' select a subset of these via the `include` argument with the following 
#' accepted values:
#' 
#' * `"effects"`: Treatment effects over time (observed series - synthetic control).
#' * `"synthetic"`: Synthetic control values over time.
#' * `"weights"`: Donor weights used in the synthetic control.
#' * `"cumulative_effects"`: Cumulative average treatment effects in the 
#'   post-treatment period.
#' * `"average_effects"`: Average treatment effects in pre- and post-treatment 
#'   periods.
#' * `"relative_change"`: Relative change from last pre-treatment synthetic 
#'    control value in the post-treatment period.
#' * `"parameters"`: Model parameters excluding the weights, e.g., 
#'   intercept (`alpha`), standard deviation of noise term (`sigma`); and in a 
#'   case of covariates, regression coefficients (`beta`).
#' * `"effective_donors"`: Effective number of donor units based on weight 
#'   concentration \eqn{1/\sum(w^2)}.
#' * `"RMSE"`: Root mean squared errors for pre- and post-treatment periods, and 
#'   their ratio.
#' * `"R2"`: Bayesian \eqn{R^2} measure of model fit.
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param include \[`character()` or `NULL`]\cr Character vector defining which 
#'   posterior summaries to include in the output. See details for accepted 
#'   values. If `NULL` (the default), all available summaries are included.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A named list of data frames containing the requested posterior 
#'   summaries.
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
  
  times <- get_times(object)
  T_pre <- get_T_pre(object)
  T_total <- get_T_total(object)
  time <- get_time(object)
  
  out <- list()
  
  if ("effects" %in% include) {
    out$effects <- treatment_effect(object, probs)
  }
  if ("synthetic" %in% include) {
    out$synthetic <- synthetic_control(object, probs)
  }
  if ("weights" %in% include) {
    out$weights <- donor_weights(object, probs)
  }
  if ("cumulative_effects" %in% include) {
    post_times <- times[(T_pre + 1L):T_total]
    out$cumulative_effects <- as_draws(object, "avg_effect_post_cumulative") |> 
      summarize_with_probs(probs) |> 
      mutate("{time}" := .env$post_times, .before = 1L) |> 
      select(-"variable")
  }
  if ("relative_change" %in% include) {
    post_times <- times[(T_pre + 1L):T_total]
    out$relative_change <- as_draws(object, "relative_change") |> 
      summarize_with_probs(probs) |> 
      mutate("{time}" := .env$post_times, .before = 1L) |> 
      select(-"variable")
  }
  if ("parameters" %in% include) {
    all_pars <- object$stanfit@model_pars
    pars <- intersect(
      c("alpha", "beta", "sigma", "tau"), all_pars
    )
    out$parameters <- as_draws(object, pars) |> 
      summarize_with_probs(probs)
  }
  if ("RMSE" %in% include) {
    out$RMSE <- rmse(object, probs)
  }
  if ("R2" %in% include) {
    out$R2 <- bayes_R2(object, probs = probs)
  }
  if ("effective_donors" %in% include) {
    out$effective_donors <- effective_donors(object, probs)
  }
  if ("average_effects" %in% include) {
    out$average_effects <- as_draws(
      object, c("avg_effect_pre", "avg_effect_post")
    ) |> 
      summarize_with_probs(probs) |> 
      mutate(
        variable = case_when(
          variable == "avg_effect_pre" ~ "Average pre-treatment effect",
          variable == "avg_effect_post" ~ "Average post-treatment effect"
        )
      )
  }
  out
}
