#' @export
#' @rdname rmse
rmse <- function(x, ...) {
  UseMethod("rmse", x)
}
#' Extract root mean squared errors of a Bayesian synthetic control model
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param average \[`logical(1)`]\cr If `TRUE`, returns the
#' average RMSEs over treated units in case of multiple treated units.
#' If `FALSE` (the default), unit-specific values are returned.
#' @param summary \[`logical(1)`]\cr If `TRUE` (the default), returns posterior
#'   mean, standard deviation, posterior quantiles (as defined by the
#'   `probs` argument), and MCMC convergence measures.
#'   If `FALSE`, returns the posterior samples instead.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior samples (`summary = FALSE`) in long format.
#' @rdname rmse
#' @aliases rmse
#' @export
#' @examples
#' rmse(fit_single_treated)
rmse.bscmfit <- function(
  x,
  average = FALSE,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  test_summary(summary)
  probs <- sort_probs(probs)
  stopifnot_(
    checkmate::test_flag(average),
    "Argument {.arg average} must be a single {.cls logical} value."
  )
  stopifnot_(
    !is.null(x$data),
    "The model fit {.arg x} does not contain the original data. You probably
    used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  time <- get_time(x)
  unit <- get_unit(x)
  N <- get_N(x)
  treatment <- get_treatment(x)
  d <- treatment_effect(x, type = "time", average = FALSE, summary = FALSE)
  cols <- c(unit, treatment)
  d <- d |>
    dplyr::summarise(
      rmse = sqrt(posterior::rvar_mean(.data$effect^2)),
      .by = dplyr::any_of(cols)
    )
  if (average && N > 1) {
    d <- d |>
      dplyr::summarise(rmse = posterior::rvar_mean(rmse), .by = .env$treatment)
  }
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(rmse, probs)) |>
      dplyr::select(-"rmse", -"variable")
  }
  if (N == 1) {
    d <- d |> dplyr::select(-any_of(unit))
  }
  d
}
