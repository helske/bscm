#' @export
#' @rdname covariate_imbalance
covariate_imbalance <- function(x, ...) {
  UseMethod("covariate_imbalance", x)
}
#' Covariate imbalance of Bayesian synthetic control model
#'
#' For models with covariates, returns the covariate imbalances
#' \deqn{\delta_{t} =
#' \sqrt{\frac{1}{K}\sum_{k=1}^K(x_{k,0,t} - \bar x_{k,0,t})^2},}
#' \eqn{t=1,\ldots,T}, where
#' \eqn{\bar x_{k,0,t} = \sum_{j=1}^J \omega_j x_{k, j, t}} and \eqn{x_{k,0,t}}
#' is the value of \eqn{k}th covariate of a treated unit at time t, and
#' similarly for donors \eqn{j=1,\ldots,J}. This is computed separately for
#' each treated unit in case of multiple treated units.
#'
#' @inheritParams bscm_postprocessing
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @rdname covariate_imbalance
#' @aliases covariate_imbalance
#' @seealso [covariate_adjustment()] for the covariate adjustments, and
#'   [plot_covariate_imbalance()] for visualizing the imbalances.
#' @export
#' @examples
#' covariate_imbalance(fit_single_treated, probs = c(0.05, 0.95))
#'
covariate_imbalance.bscmfit <- function(
  x,
  average = FALSE,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  warn_deprecated_plot(..., replacement = "plot_covariate_imbalance")
  check_flag(average, "average")
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  check_has_predictors(x)
  d <- covariate_imbalance_draws(x, average)
  if (summary) {
    d <- summarise_column(d, "imbalance", probs)
  }
  d
}

#' Posterior draws of the covariate imbalances of a BSCM
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param average \[`logical(1)`]\cr If `TRUE`, the imbalances are averaged
#'   over the treated units within each posterior draw.
#' @noRd
covariate_imbalance_draws <- function(x, average = FALSE) {
  X <- get_Xs(x)
  time <- get_time(x)
  times <- get_times(x)
  d <- lapply(
    seq_len(get_N(x)),
    \(i) {
      dplyr::tibble(
        "{time}" := times,
        imbalance = sqrt(
          posterior::rvar_apply(
            covariate_delta(x, i, X)^2,
            1L,
            posterior::rvar_mean
          )
        )
      )
    }
  ) |>
    stats::setNames(get_treated(x)) |>
    dplyr::bind_rows(.id = get_unit(x))
  if (average && get_N(x) > 1L) {
    d <- d |>
      dplyr::summarise(
        imbalance = posterior::rvar_mean(.data$imbalance),
        .by = dplyr::all_of(time)
      )
  }
  d
}
