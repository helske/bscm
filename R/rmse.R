#' @export
#' @rdname rmse
rmse <- function(x, ...) {
  UseMethod("rmse", x)
}
#' Extract root mean squared errors of a Bayesian synthetic control model
#'
#' Returns posterior draws (or summaries) of root mean squared errors (RMSEs)
#' for the pre-treatment and post-treatment periods.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param average \[`logical(1)`]\cr If `TRUE`, returns the
#' average RMSEs over treated units in case of multiple treated units.
#' If `FALSE` (the default), unit-specific values are returned.
#' @param summary \[`logical(1)`]\cr If `TRUE` (the default), returns posterior
#'   mean, standard deviation, posterior quantiles (as defined by the
#'   `probs` argument), and MCMC convergence measures.
#'   If `FALSE`, returns the posterior draws instead.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`).
#' @rdname rmse
#' @aliases rmse
#' @export
#' @examples
#' rmse(fit_single_treated, probs = c(0.01, 0.5, 0.8))
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
  d <- treatment_effect(x, average = FALSE, summary = FALSE)
  cols <- c(unit, treatment)
  d <- d |>
    dplyr::summarise(
      rmse = sqrt(posterior::rvar_mean(.data$effect^2)),
      .by = dplyr::any_of(cols)
    )
  if (average && N > 1) {
    d <- d |>
      dplyr::summarise(
        rmse = posterior::rvar_mean(.data$rmse),
        .by = .env$treatment
      )
  }
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(.data$rmse, probs)) |>
      dplyr::select(-"rmse", -"variable")
  }
  if (N == 1) {
    d <- d |> dplyr::select(-any_of(unit))
  }
  d
}
#' @export
#' @rdname rmse_ratio
rmse_ratio <- function(x, ...) {
  UseMethod("rmse_ratio", x)
}

#' Extract the ratio of post- to pre-treatment RMSE
#'
#' Returns posterior draws (or summaries) of the ratio of the
#' post-treatment RMSE to the pre-treatment RMSE.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param average \[`logical(1)`]\cr If `TRUE`, returns the
#' average ratio over treated units in case of multiple treated units.
#' If `FALSE` (the default), unit-specific values are returned.
#' @param summary \[`logical(1)`]\cr If `TRUE` (the default), returns posterior
#' summaries. If `FALSE`, returns posterior draws.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#' Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A tibble of posterior summaries (`summary = TRUE`) or
#' posterior draws (`summary = FALSE`).
#' @export
#' @rdname rmse_ratio
#' @examples
#' rmse_ratio(fit_single_treated, probs = c(0.01, 0.5, 0.8))
rmse_ratio.bscmfit <- function(
  x,
  average = FALSE,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  test_summary(summary)
  probs <- sort_probs(probs)
  unit <- get_unit(x)
  treatment <- get_treatment(x)
  d <- rmse(
    x,
    average = average,
    summary = FALSE
  )

  d <- d |>
    dplyr::summarise(
      ratio = .data$rmse[.data[[treatment]] == 1] /
        .data$rmse[.data[[treatment]] == 0],
      .by = dplyr::any_of(unit)
    )

  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(.data$ratio, probs)) |>
      dplyr::select(-"ratio", -"variable")
  }

  if (get_N(x) == 1) {
    d <- d |> dplyr::select(-dplyr::any_of(unit))
  }
  d
}
