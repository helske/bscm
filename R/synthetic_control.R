#' @export
#' @rdname synthetic_control
synthetic_control <- function(x, ...) {
  UseMethod("synthetic_control", x)
}
#' Extract synthetic control series of a Bayesian synthetic control model
#'
#' @inheritParams rmse.bscmfit
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or
#'   posterior samples (`summary = FALSE`) in long format.
#' @rdname synthetic_control
#' @aliases synthetic_control
#' @export
#' @examples
#' synthetic_control(fit_single_treated) |> tail()
synthetic_control.bscmfit <- function(
  x,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  test_summary(summary)
  probs <- sort_probs(probs)
  time <- get_time(x)
  times <- get_times(x)
  N <- get_N(x)
  for_plots <- list(...)$for_plots %||% FALSE
  unit <- get_unit(x)
  treated <- get_treated(x)
  values <- as_draws(x, "y_rep")
  units <- rep(treated, each = length(times))
  times2 <- rep(times, times = length(treated))
  format_posterior_output(
    values,
    summary = summary,
    probs = probs,
    variable = "Synthetic control",
    for_plots = for_plots
  ) |>
    add_output_column(
      name = time,
      values = times2,
      summary = summary
    ) |>
    add_output_column(
      name = unit,
      values = units,
      summary = summary
    )
}
