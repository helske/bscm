#' @export
#' @rdname donor_weights
donor_weights <- function(x, ...) {
  UseMethod("donor_weights", x)
}
#' Extract donor weights of a Bayesian synthetic control model
#'
#' @inheritParams rmse.bscmfit
#' @param ... Ignored.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior samples (`summary = FALSE`) in long format.
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
  test_summary(summary)
  probs <- sort_probs(probs)
  unit <- get_unit(x)
  treated <- get_treated(x)
  N <- get_N(x)
  donors <- get_donors(x)
  J <- length(donors)
  d <- dplyr::tibble(
    treated = rep(treated, times = J),
    donor = rep(donors, each = N),
    weight = c(posterior::as_draws_rvars(as_draws(x, "omega"))$omega)
  )
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(.data$weight, probs)) |>
      dplyr::select(-"weight", -"variable")
  }
  d
}
