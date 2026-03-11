#' @export
#' @rdname synthetic_control
synthetic_control <- function(x, ...) {
  UseMethod("synthetic_control", x)
}
#' Extract synthetic control series of a Bayesian synthetic control model
#'
#' @export
#' @rdname synthetic_control
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.5, 0.975)`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the synthetic_control series.
synthetic_control.bscmfit <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
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
  as_draws(x, "synthetic_y") |> 
  summarize_with_probs(probs) |> 
    mutate("{time}" := .env$times, .before = 1L) |> 
    select(-"variable")
}
