#' Extract residual standard deviations of a Bayesian synthetic control
#' model
#'
#' @inheritParams rmse.bscmfit
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or 
#'   posterior samples (`summary = FALSE`) in long format.
#' @aliases sigma
#' @export
sigma.bscmfit <- function(object, summary = TRUE, 
                          probs = c(0.025, 0.975), ...) {
  test_summary(summary)
  probs <- sort_probs(probs)
  treated <- get_treated(object)
  unit <- get_unit(object)
  format_posterior_output(
    as_draws(object, "sigma"),
    summary = summary,
    probs = probs,
    variable = "sigma"
  ) |>
    add_output_column(
      name = unit,
      values = treated,
      summary = summary
    )
}
