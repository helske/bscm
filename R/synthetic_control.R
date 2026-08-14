#' @export
#' @rdname synthetic_control
synthetic_control <- function(x, ...) {
  UseMethod("synthetic_control", x)
}
#' Synthetic control series of a Bayesian synthetic control model
#'
#' Returns posterior draws (or summaries) of the pre- and post-treatment
#' trajectories of treated units, i.e. draws `y_rep`from the posterior
#' predictive distribution of the model.
#'
#' @inheritParams rmse.bscmfit
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
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
  probs <- sort_probs(probs)
  test_summary(summary)
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
  unit <- get_unit(x)
  treated <- get_treated(x)
  treatment <- get_treatment(x)
  treatments <- unlist(
    lapply(treated, \(i) rep(0:1, times = c(T_pre[i], T_total - T_pre[i])))
  )
  y_rep <- posterior_predict(x)
  d <- dplyr::tibble(
    "{unit}" := rep(treated, each = T_total),
    "{time}" := rep(times, times = N),
    "{treatment}" := treatments,
    y_rep = c(posterior::as_draws_rvars(y_rep)$y_rep)
  )
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(.data$y_rep, probs, for_plots)) |>
      dplyr::select(-"y_rep", -"variable")
  }
  if (N == 1) {
    d <- d |> dplyr::select(-dplyr::all_of(unit))
  }
  d
}
