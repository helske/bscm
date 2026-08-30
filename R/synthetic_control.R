#' @export
#' @rdname synthetic_control
synthetic_control <- function(x, ...) {
  UseMethod("synthetic_control", x)
}
#' Synthetic control series of a Bayesian synthetic control model
#'
#' Returns posterior draws (or summaries) of the pre- and post-treatment
#' trajectories of treated units, i.e. draws `y_rep` from the posterior
#' predictive distribution of the model.
#'
#' @inheritParams bscm_postprocessing
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @rdname synthetic_control
#' @aliases synthetic_control
#' @seealso [bscm::posterior_predict.bscmfit()] for the `y_rep` draws in matrix format.
#' @export
#' @examples
#' synthetic_control(fit_single_treated) |> tail()
synthetic_control.bscmfit <- function(
  x,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  check_has_data(x, "x")
  d <- sc_draws(x)
  if (summary) {
    d <- summarise_column(d, "y_rep", probs)
  }
  d
}
#' Posterior draws of the synthetic control series
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @noRd
sc_draws <- function(x) {
  y_rep <- gq_to_rvar(x, posterior_predict(x))
  treated_grid(x) |>
    dplyr::mutate(y_rep = c(.env$y_rep))
}
