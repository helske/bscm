#' @export
#' @rdname donor_weights
donor_weights <- function(x, ...) {
  UseMethod("donor_weights", x)
}
#' Extract donor weights of a Bayesian synthetic control model
#'
#' @inheritParams bscm_postprocessing
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @rdname donor_weights
#' @aliases donor_weights
#' @export
#' @examples
#' donor_weights(fit_single_treated) |> head(5)
donor_weights.bscmfit <- function(
  x,
  summary = TRUE,
  probs = c(0.025, 0.5, 0.975),
  ...
) {
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  d <- weights_draws(x)
  if (summary) {
    d <- summarise_column(d, "weight", probs)
  }
  d
}

#' Posterior draws of the donor weights of a BSCM
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @noRd
weights_draws <- function(x) {
  donors <- get_donors(x)
  J <- length(donors)
  dplyr::tibble(
    treated = rep(get_treated(x), each = J),
    donor = rep(donors, times = get_N(x)),
    weight = c(rvars_of(x, "omega"))
  )
}
