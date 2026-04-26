#' @export
#' @rdname synthetic_control
synthetic_control <- function(x, ...) {
  UseMethod("synthetic_control", x)
}
#' Extract synthetic control series of a Bayesian synthetic control model
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of synthetic controls.
#' @rdname synthetic_control
#' @aliases synthetic_control
#' @export
synthetic_control.bscmfit <- function(x, probs = c(0.025, 0.975), ...) {
  probs <- sort_probs(probs)
  time <- get_time(x)
  times <- get_times(x)
  N <- get_N(x)
  for_plots <- list(...)$for_plots %||% FALSE
  unit <- get_unit(x)
  treated <- get_treated(x)
  d <- as_draws(x, "y_rep")
  d |>
    summarise_with_probs(probs, for_plots) |>
    mutate(
      "{unit}" := rep(treated, each = length(times)),
      "{time}" := rep(times, times = length(treated)),
      .before = 1L
    ) |>
    mutate(variable = "Synthetic control")
}
