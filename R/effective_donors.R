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
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param average \[`logical(1)`]\cr If `TRUE` (the default) and the model 
#'   contains multiple treated units, an additional row with the posterior 
#'   summary of the average number of effective donors across treated units is
#'   appended to the output.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the estimated effective 
#' donors.
#' @rdname effective_donors
#' @aliases effective_donors
#' @export
effective_donors.bscmfit <- function(x, probs = c(0.025, 0.975),
                                     average = TRUE, ...) {
  
  test_probs(probs)
  stopifnot_(
    checkmate::test_flag(average),
    "Argument {.arg average} must be a single {.cls logical} value."
  )
  N <- get_N(x)
  if (average && N > 1) {
    as_draws(x, "avg_effective_donors") |> 
      summarise_with_probs(probs) |> 
      mutate(variable = "Average effective donors")
  } else {
    unit <- get_unit(x)
    treated <- get_treated(x)
    as_draws(x, "effective_donors") |> 
      summarise_with_probs(probs) |>
      mutate("{unit}" := treated, .before = 1L) |> 
      mutate(variable = "Effective donors")
  }
}
