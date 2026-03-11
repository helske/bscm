#' @export
#' @rdname treatment_effect
treatment_effect <- function(x, ...) {
  UseMethod("treatment_effect", x)
}
#' Treatment effect estimates of a Bayesian synthetic control model
#'
#' @export
#' @rdname treatment_effect
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the treatment effects per 
#' time point.
treatment_effect.bscmfit <- function(x, probs = c(0.025, 0.975), ...) {
  
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
  time <- get_time(x)
  times <- get_times(x)
  as_draws(x, "effect") |> 
  summarize_with_probs(probs) |> 
    mutate("{time}" := .env$times, .before = 1L) |> 
    select(-"variable")
}
