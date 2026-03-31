#' @export
#' @rdname treatment_effect
treatment_effect <- function(x, ...) {
  UseMethod("treatment_effect", x)
}
#' Treatment effect estimates of a Bayesian synthetic control model
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param average \[`logical(1)`]\cr If `TRUE` (the default), returns the 
#' average treatment effects over treated units (per time point) in case of 
#' multiple treated units. If `FALSE`, unit-specific effects are returned.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the treatment effects per 
#' time point.
#' @rdname treatment_effect
#' @aliases treatment_effect
#' @export
treatment_effect.bscmfit <- function(x, probs = c(0.025, 0.975), 
                                     average = TRUE, ...) {
  
  test_probs(probs)
  stopifnot_(
    checkmate::test_flag(average),
    "Argument {.arg average} must be a single {.cls logical} value."
  )
  time <- get_time(x)
  times <- get_times(x)
  N <- get_N(x)
  for_plots <- list(...)$for_plots %||% FALSE
  d <- as_draws(x, "effect")
  
  if (average && N > 1) {
    out <- as_draws_rvars(d)[[1]] |> 
      rvar_apply(1, rvar_mean) |> 
      summarise_with_probs(probs, for_plots) |>
      mutate(
        "{time}" := times,
        .before = 1L
      ) |>
      mutate(variable = "Average treatment effect")
  } else {
    unit <- get_unit(x)
    treated <- get_treated(x)
    out <- d |>
      summarise_with_probs(probs, for_plots) |>
      mutate(
        "{unit}" := rep(treated, each = length(times)),
        "{time}" := rep(times, times = length(treated)),
        .before = 1L
      ) |> 
      mutate(variable = "Treatment effect")
  }
  out
}
