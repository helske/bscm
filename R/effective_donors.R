#' @export
#' @rdname effective_donors
effective_donors <- function(x, ...) {
  UseMethod("effective_donors", x)
}
#' Extract the number of effective donors in a Bayesian synthetic control model
#'
#' Effective donors is defined as \eqn{1 / \sum_{j=1}^J \omega_j^2},
#' where \eqn{\omega_j} is the donor weight of control unit \eqn{j}.
#' @export
#' @rdname effective_donors
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.5, 0.975)`.
#' @return A `data.frame` of posterior summaries of the estimated effective 
#' donors.
effective_donors.bscmfit <- function(x, probs = c(0.025, 0.975), ...) {
 
  stopifnot_(
    checkmate::test_numeric(
      probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values between
     0 and 1."
  )
  as_draws(x, "effective_donors") |> 
    summarize_with_probs(probs)
}
