#' @export
#' @rdname donor_weights
donor_weights <- function(x, ...) {
  UseMethod("donor_weights", x)
}
#' Extract donor weights of a Bayesian synthetic control model
#'
#' @inheritParams rmse.bscmfit
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the donor weights
#' @rdname donor_weights
#' @aliases donor_weights
#' @export
donor_weights.bscmfit <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
  probs <- sort_probs(probs)
  donors <- get_donors(x)
  unit <- get_unit(x)
  treated <- get_treated(x)
  N <- get_N(x)
  J <- get_J(x)
  out <- lapply(
    seq_len(N),
    \(i) {
      pars <- paste0("omega[", i, ",", seq_len(J), "]")
      as_draws(x, pars) |>
        summarise_with_probs(probs) |>
        mutate(
          treated_unit = .env$treated[i],
          "{unit}" := .env$donors,
          .before = 1L
        ) |>
        select(-"variable")
    }
  )
  bind_rows(out)
}
