#' @export
#' @rdname donor_weights
donor_weights <- function(x, ...) {
  UseMethod("donor_weights", x)
}
#' Extract donor weights of a Bayesian synthetic control model
#'
#' @inheritParams rmse.bscmfit
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or
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
  donors <- get_donors(x)
  unit <- get_unit(x)
  treated <- get_treated(x)
  N <- get_N(x)
  J <- get_J(x)
  out <- lapply(
    seq_len(N),
    \(i) {
      pars <- paste0("omega[", i, ",", seq_len(J), "]")
      format_posterior_output(
        as_draws(x, pars),
        summary = summary,
        probs = probs,
        variable = "Donor weight"
      ) |>
        add_output_column(
          name = unit,
          values = donors,
          summary = summary
        ) |>
        add_output_column(
          name = "treated_unit",
          values = treated[i],
          summary = summary
        )
    }
  )
  bind_rows(out)
}
