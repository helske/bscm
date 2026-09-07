#' Summarise posterior draws of an estimated Bayesian synthetic control model
#'
#' Returns posterior summaries of the model parameters, i.e. the posterior
#' mean, standard deviation, quantiles as defined by the `probs` argument, and
#' the MCMC convergence measures. The parameters are named as in the posterior
#' draws, see [as_draws.bscmfit()].
#'
#' Note that many (more) summaries are available via separate methods, e.g.,
#' [donor_weights()], [bayes_R2()], [rmse()], and
#' [average_treatment_effect()], and that printing the model fit shows the key
#' quantities of the model, see [print.bscmfit()].
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param parameters \[`character()`]\cr Vector of parameter names. When
#'   `NULL` (the default), corresponds to a relevant subset of the model
#'   parameters, as in [as_draws.bscmfit()].
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#'   Default is `c(0.025, 0.975)`.
#' @param diagnostics \[`logical(1)`]\cr If `TRUE` (the default), the output
#'   also contains the MCMC convergence measures and the Monte Carlo standard
#'   error of the mean.
#' @param ... Additional arguments passed to [as_draws.bscmfit()].
#' @return A `tibble` with posterior summaries of the model parameters.
#' @seealso [as_draws.bscmfit()] for the posterior draws themselves.
#' @aliases summary
#' @export
#' @examples
#' summary(fit_single_treated)
summary.bscmfit <- function(
  object,
  parameters = NULL,
  probs = c(0.025, 0.975),
  diagnostics = TRUE,
  ...
) {
  probs <- sort_probs(probs)
  check_flag(diagnostics, "diagnostics")
  object |>
    as_draws(parameters = parameters, ...) |>
    summarise_with_probs(probs, diagnostics)
}
