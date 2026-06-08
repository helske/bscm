#' @export
#' @rdname effective_donors
effective_donors <- function(x, ...) {
  UseMethod("effective_donors", x)
}
#' Extract the number of effective donors in a Bayesian synthetic control model
#'
#' Effective donors is defined as \eqn{1 / \sum_{j=1}^J \omega_j^2},
#' where \eqn{\omega_j} is the donor weight of control unit \eqn{j}.
#'
#' @inheritParams rmse.bscmfit
#' @param average \[`logical(1)`]\cr If `TRUE`, returns the
#' average effective donors over treated units in case of
#' multiple treated units. The default is `FALSE`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or
#'   posterior samples (`summary = FALSE`) in long format.
#' @rdname effective_donors
#' @aliases effective_donors
#' @export
#' @examples
#' effective_donors(fit_single_treated)
effective_donors.bscmfit <- function(
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

  treated <- get_treated(x)
  unit <- get_unit(x)
  N <- get_N(x)
  omega <- posterior::as_draws_rvars(as_draws(x, "omega"))$omega
  d <- dplyr::tibble(
    "{unit}" := treated,
    ess = posterior::rvar_apply(omega, 1, \(x) 1 / posterior::rvar_sum(x^2))
  )
  if (average && N > 1) {
    d <- d |>
      dplyr::summarise(ess = posterior::rvar_mean(.data$ess))
  }
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(.data$ess, probs)) |>
      dplyr::select(-"ess", -"variable")
  }
  if (N == 1) {
    d <- d |> dplyr::select(-dplyr::all_of(unit))
  }
  d
}
