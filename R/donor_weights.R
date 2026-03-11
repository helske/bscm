#' @export
#' @rdname donor_weights
donor_weights <- function(x, ...) {
  UseMethod("donor_weights", x)
}
#' Extract donor weights of a Bayesian synthetic control model
#'
#' @export
#' @rdname donor_weights
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.5, 0.975)`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the donor weights
donor_weights.bscmfit <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
 
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
  donors <- get_donors(x)
  unit <- get_unit(x)
  as_draws(x, "omega") |> 
    summarize_with_probs(probs) |> 
    mutate("{unit}" := .env$donors, .before = 1L) |> 
    select(-"variable")
}
