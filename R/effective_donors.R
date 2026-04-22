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
#' @param average \[`logical(1)`]\cr If `TRUE` (the default), returns the 
#' average effective donors over treated units in case of 
#' multiple treated units.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the estimated effective 
#' donors.
#' @rdname effective_donors
#' @aliases effective_donors
#' @export
effective_donors.bscmfit <- function(x, average = TRUE, probs = c(0.025, 0.975),
                                      ...) {
  
  test_probs(probs)
  stopifnot_(
    checkmate::test_flag(average),
    "Argument {.arg average} must be a single {.cls logical} value."
  )
  
  N <- get_N(x)
  omega <- as_draws_rvars(as_draws(x, "omega"))$omega
  eff <- rvar_apply(omega,1,\(x) 1/rvar_sum(x^2))
 
  if (average && N > 1) {
    rvar_mean(eff) |>
      summarise_with_probs(probs) |>
      mutate(variable = "Average effective donors")
  } else {
    unit <- get_unit(x)
    treated <- get_treated(x)
    eff |>
      summarise_with_probs(probs) |>
      mutate("{unit}" := treated, .before = 1L) |>
      mutate(variable = "Effective donors")
  }
}
