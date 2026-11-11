#' Summarise posterior draws of an estimated Bayesian synthetic control model
#'
#' Returns posterior summaries of the model parameters, i.e. the posterior
#' mean, standard deviation, quantiles as defined by the `probs` argument, and
#' the MCMC convergence measures, see [bscm_prior()] for the parameter
#' names.
#'
#' Note that many (more) summaries are available via separate methods, e.g.,
#' [donor_weights()], [bayes_R2()], [rmse()], and
#' [average_treatment_effect()], and that printing the model fit shows the key
#' quantities of the model, see [print.bscmfit()].
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param parameters \[`character()`]\cr Names of the model parameters to
#'   summarise, among the parameters of the model listed in [bscm_prior()] other
#'   than weights, which are summarised by [donor_weights()]. When `NULL` (the
#'   default), all of these parameters. The other quantities of the model, such
#'   as the time-varying parts of the coefficients and the spline coefficients,
#'   are available via [as_draws.bscmfit()].
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#'   Default is `c(0.025, 0.975)`.
#' @param diagnostics \[`logical(1)`]\cr If `TRUE` (the default), the output
#'   also contains the MCMC convergence measures and the Monte Carlo standard
#'   error of the mean.
#' @param ... Ignored.
#' @return A `tibble` with one row per element of each model parameter, with
#'   columns `parameter`, `term` (predictor or tilt direction), a
#'   column named by the unit variable (treated unit), and the posterior
#'   summaries. Columns `term` and the unit column are omitted when not
#'   needed.
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
  available <- c(
    setdiff(names(object$priors), "weights"),
    if (has_zeta(object)) "concentration",
    if (!is.null(object$priors$weights$tilt)) "tilt"
  )
  parameters <- parameters %||% available
  stopifnot_(
    checkmate::test_subset(parameters, available, empty.ok = FALSE),
    c(
      "Argument {.arg parameters} must contain parameters of the model.",
      i = "Parameters of this model are {.val {available}}."
    )
  )
  parameter_draws(object, intersect(available, parameters)) |>
    summarise_column("value", probs, diagnostics)
}
