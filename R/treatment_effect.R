#' @export
#' @rdname treatment_effect
treatment_effect <- function(x, ...) {
  UseMethod("treatment_effect", x)
}
#' Treatment effect estimates of a Bayesian synthetic control model
#'
#' @inheritParams bscm_postprocessing
#' @param average \[`logical(1)`]\cr If `TRUE`, returns the average treatment
#'   effect over treated units at each time since treatment (event time). If
#'   `FALSE` (the default), unit-specific effects at each calendar time point
#'   are returned.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @rdname treatment_effect
#' @aliases treatment_effect
#' @export
#' @examples
#' treatment_effect(fit_single_treated) |> tail()
treatment_effect.bscmfit <- function(
  x,
  average = FALSE,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  check_flag(average, "average")
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  check_has_data(x, "x")
  d <- effects_draws(x, average = average && get_N(x) > 1L)
  if (summary) {
    d <- summarise_column(d, "effect", probs)
  }
  d
}
#' @export
#' @rdname average_treatment_effect
average_treatment_effect <- function(x, ...) {
  UseMethod("average_treatment_effect", x)
}
#' Average treatment effects for pre- and post-treatment periods of a BSCM
#'
#' @inheritParams treatment_effect.bscmfit
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @rdname average_treatment_effect
#' @aliases average_treatment_effect
#' @export
#' @examples
#' average_treatment_effect(fit_single_treated)
average_treatment_effect.bscmfit <- function(
  x,
  average = FALSE,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  check_flag(average, "average")
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  check_has_data(x, "x")
  N <- get_N(x)
  unit <- get_unit(x)
  treatment <- get_treatment(x)
  average <- average && N > 1
  cols <- c(if (!average) unit, treatment)
  d <- effects_draws(x) |>
    dplyr::summarise(
      effect = posterior::rvar_mean(.data$effect),
      .by = dplyr::all_of(cols)
    )
  if (summary) {
    d <- summarise_column(d, "effect", probs)
  }
  d
}
#' Posterior draws of the treatment effects of a BSCM
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param average \[`logical(1)`]\cr If `TRUE`, the effects are averaged over
#'   the treated units within each posterior draw at each time point since
#'   treatment.
#' @noRd
effects_draws <- function(x, average = FALSE) {
  y_rep <- gq_to_rvar(x, posterior_predict(x))
  effect <- get_stan_y(x) - y_rep
  d <- treated_grid(x, event_time = TRUE) |>
    dplyr::mutate(effect = c(.env$effect))
  if (average) {
    d <- d |>
      dplyr::summarise(
        effect = posterior::rvar_mean(.data$effect),
        .by = "time_since_treatment"
      ) |>
      dplyr::arrange(.data$time_since_treatment) |>
      dplyr::rename("{get_time(x)}" := "time_since_treatment")
  } else {
    d <- dplyr::select(d, -"time_since_treatment")
  }
  d
}
