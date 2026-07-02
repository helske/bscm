#' @noRd
effects_draws <- function(x) {
  time <- get_time(x)
  times <- get_times(x)
  T_pre <- get_T_pre(x)
  T_total <- get_T_total(x)
  N <- get_N(x)
  unit <- get_unit(x)
  treated <- get_treated(x)
  treatment <- get_treatment(x)
  y_rep <- posterior_predict(x)
  effect <- get_stan_y(x) - posterior::as_draws_rvars(y_rep)$y_rep
  treatments <- unlist(
    lapply(treated, \(i) rep(0:1, times = c(T_pre[i], T_total - T_pre[i])))
  )
  event_times <- unlist(
    lapply(treated, \(i) seq_len(T_total) - T_pre[i] - 1L)
  )
  tibble(
    "{unit}" := rep(treated, each = T_total),
    "{time}" := rep(times, times = N),
    time_since_treatment = event_times,
    "{treatment}" := treatments,
    effect = c(effect)
  )
}

#' @export
#' @rdname treatment_effect
treatment_effect <- function(x, ...) {
  UseMethod("treatment_effect", x)
}
#' Treatment effect estimates of a Bayesian synthetic control model
#'
#' @inheritParams rmse.bscmfit
#' @param average \[`logical(1)`]\cr If `TRUE` (the default), returns the
#'   average treatment effect over treated units at each time since treatment
#'   (event time). If `FALSE`, unit-specific effects at each calendar time
#'   point are returned.
#' @param ... Ignored.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @rdname treatment_effect
#' @aliases treatment_effect
#' @export
#' @examples
#' treatment_effect(fit_single_treated) |> tail()
treatment_effect.bscmfit <- function(
    x,
    average = TRUE,
    summary = TRUE,
    probs = c(0.025, 0.975),
    ...
) {
  probs <- sort_probs(probs)
  test_summary(summary)
  stopifnot_(
    checkmate::test_flag(average),
    "Argument {.arg average} must be a single {.cls logical} value."
  )
  stopifnot_(
    !is.null(x$data),
    "The model fit {.arg x} does not contain the original data. You probably
    used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  for_plots <- list(...)$for_plots %||% FALSE
  N <- get_N(x)
  unit <- get_unit(x)
  treatment <- get_treatment(x)
  average <- average && N > 1
  d <- effects_draws(x)
  if (average) {
    d <- d |>
      dplyr::summarise(
        effect = posterior::rvar_mean(effect),
        .by = "time_since_treatment"
      ) |>
      dplyr::arrange(time_since_treatment)
  } else {
    d <- dplyr::select(d, -"time_since_treatment")
  }
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(effect, probs, for_plots)) |>
      dplyr::select(-"effect", -"variable")
  }
  if (N == 1) {
    d <- d |> dplyr::select(-dplyr::any_of(unit))
  }
  d
}

#' @export
#' @rdname average_treatment_effect
average_treatment_effect <- function(x, ...) {
  UseMethod("average_treatment_effect", x)
}
#' Pre- and post-treatment average effects of a Bayesian synthetic control
#' model
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
    average = TRUE,
    summary = TRUE,
    probs = c(0.025, 0.975),
    ...
) {
  probs <- sort_probs(probs)
  test_summary(summary)
  stopifnot_(
    checkmate::test_flag(average),
    "Argument {.arg average} must be a single {.cls logical} value."
  )
  stopifnot_(
    !is.null(x$data),
    "The model fit {.arg x} does not contain the original data. You probably
    used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  N <- get_N(x)
  unit <- get_unit(x)
  treatment <- get_treatment(x)
  average <- average && N > 1
  cols <- c(if (!average) unit, treatment)
  d <- effects_draws(x) |>
    dplyr::summarise(
      effect = posterior::rvar_mean(effect),
      .by = dplyr::all_of(cols)
    )
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(effect, probs)) |>
      dplyr::select(-"effect", -"variable")
  }
  if (N == 1) {
    d <- d |> dplyr::select(-dplyr::any_of(unit))
  }
  d
}
