#' @export
#' @rdname treatment_effect
treatment_effect <- function(x, ...) {
  UseMethod("treatment_effect", x)
}
#' Treatment effect estimates of a Bayesian synthetic control model
#'
#' @inheritParams rmse.bscmfit
#' @param type \[`character(1)`]\cr Type of treatment effect to compute.
#'   `"time"` (the default) returns effects at each time point.
#'   `"average"` returns the average pre- and post-treatment effects.
#' @param average \[`logical(1)`]\cr If `TRUE` (the default), returns the
#' average effects over treated units in case of
#' multiple treated units. If `FALSE`, unit-specific effects are returned.
#' @param ... Ignored.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior samples (`summary = FALSE`) in long format.
#' @rdname treatment_effect
#' @aliases treatment_effect
#' @export
#' @examples
#' treatment_effect(fit_single_treated) |> tail()
treatment_effect.bscmfit <- function(
  x,
  type = "time",
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
  type <- try_(match.arg(type, c("time", "average")))
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val time} or {.val average}."
  )
  stopifnot_(
    !is.null(x$data),
    "The model fit {.arg x} does not contain the original data. You probably
    used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  for_plots <- list(...)$for_plots %||% FALSE
  time <- get_time(x)
  times <- get_times(x)
  T_pre <- get_T_pre(x)
  T_total <- get_T_total(x)
  N <- get_N(x)
  average <- average && N > 1
  unit <- get_unit(x)
  treated <- get_treated(x)
  treatment <- get_treatment(x)
  y_rep <- posterior_predict(x)
  effect <- get_stan_y(x) - posterior::as_draws_rvars(y_rep)$y_rep
  treatments <- unlist(
    lapply(treated, \(i) rep(0:1, times = c(T_pre[i], T_total - T_pre[i])))
  )
  d <- tibble(
    "{unit}" := rep(treated, each = T_total),
    "{time}" := rep(times, times = N),
    "{treatment}" := treatments,
    effect = c(effect)
  )
  if (type == "time" && average) {
    d <- d |>
      dplyr::summarise(
        effect = posterior::rvar_mean(effect),
        .by = dplyr::all_of(c(time, treatment))
      )
  }
  if (type == "average") {
    cols <- c(if (!average) unit, treatment)
    d <- d |>
      dplyr::summarise(
        effect = posterior::rvar_mean(effect),
        .by = dplyr::all_of(cols)
      )
  }

  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(effect, probs, for_plots)) |>
      dplyr::select(-"effect", -"variable")
  }
  if (N == 1) {
    d <- d |> dplyr::select(-any_of(unit))
  }
  d
}
