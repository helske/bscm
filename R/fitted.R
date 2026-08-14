#' Expected values of posterior predictive distribution of Bayesian synthetic
#' control model
#'
#' Returns posterior draws (or summaries) of expected values `y_mean` of the
#' posterior predictive distribution of the model.
#'
#' @inheritParams rmse.bscmfit
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @aliases fitted
#' @export
#' @examples
#' fitted(fit_single_treated) |> head()
fitted.bscmfit <- function(
  object,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  probs <- sort_probs(probs)
  test_summary(summary)
  stopifnot_(
    !is.null(object$data),
    "The model fit {.arg object} does not contain the original data. You
    probably used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  for_plots <- list(...)$for_plots %||% FALSE
  time <- get_time(object)
  times <- get_times(object)
  T_pre <- get_T_pre(object)
  T_total <- get_T_total(object)
  N <- get_N(object)
  unit <- get_unit(object)
  treated <- get_treated(object)
  treatment <- get_treatment(object)
  treatments <- unlist(
    lapply(treated, \(i) rep(0:1, times = c(T_pre[i], T_total - T_pre[i])))
  )
  d <- dplyr::tibble(
    "{unit}" := rep(treated, each = T_total),
    "{time}" := rep(times, times = N),
    "{treatment}" := treatments,
    y_mean = c(posterior::as_draws_rvars(object$y_mean)$y_mean)
  )
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(.data$y_mean, probs, for_plots)) |>
      dplyr::select(-"y_mean", -"variable")
  }
  if (N == 1) {
    d <- d |> dplyr::select(-dplyr::all_of(unit))
  }
  d
}
