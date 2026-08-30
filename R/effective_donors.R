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
#' @inheritParams bscm_postprocessing
#' @param average \[`logical(1)`]\cr If `TRUE`, returns the
#' average effective donors over treated units in case of
#' multiple treated units. The default is `FALSE`.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
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
  check_flag(average, "average")
  check_flag(summary, "summary")
  probs <- sort_probs(probs)

  unit <- get_unit(x)
  omega <- rvars_of(x, "omega")
  d <- dplyr::tibble(
    "{unit}" := get_treated(x),
    ess = posterior::rvar_apply(omega, 2, \(x) 1 / posterior::rvar_sum(x^2))
  )
  if (average && get_N(x) > 1L) {
    d <- d |>
      dplyr::summarise(ess = posterior::rvar_mean(.data$ess))
  }
  if (summary) {
    d <- summarise_column(d, "ess", probs)
  }
  d
}
