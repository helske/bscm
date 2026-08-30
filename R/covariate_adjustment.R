#' @export
#' @rdname covariate_adjustment
covariate_adjustment <- function(x, ...) {
  UseMethod("covariate_adjustment", x)
}
#' Covariate adjustments for Bayesian synthetic control model
#'
#' Computes the posterior distribution of the covariate adjustment to the
#' synthetic control for each predictor, defined as the
#' regression coefficient multiplied by the difference between the predictor
#' value of the treated unit and the corresponding weighted predictor value
#' of the donor units.
#'
#' For predictor \eqn{k}, treated unit \eqn{i}, and time point \eqn{t}, the
#' covariate adjustment is
#' \deqn{
#'   b_{tk}
#'   \left(
#'     X_{y,itk} - \sum_j \omega_{ji} X_{z,jtk}
#'   \right),
#' }
#' where \eqn{b_{tk} = \beta_k} for predictors with a fixed coefficient and
#' \eqn{b_{tk} = \beta_k + \gamma_{tk}} for predictors with a time-varying
#' coefficient.
#'
#' If there are multiple treated units, the covariate adjustments are
#' returned separately for each treated unit, unless `average = TRUE`.
#'
#' @inheritParams bscm_postprocessing
#' @param x \[`bscmfit`]\cr An object of class `bscmfit`.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) of the covariate adjustment for each
#'   predictor and time point.
#' @rdname covariate_adjustment
#' @aliases covariate_adjustment
#' @seealso [covariate_imbalance()] for the covariate imbalances, and
#'   [plot_covariate_adjustment()] for visualizing the adjustments.
#' @export
#' @examples
#' covariate_adjustment(fit_single_treated) |> head()
#'
covariate_adjustment.bscmfit <- function(
  x,
  average = FALSE,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  warn_deprecated_plot(..., replacement = "plot_covariate_adjustment")
  check_flag(average, "average")
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  check_has_predictors(x)
  d <- covariate_adjustment_draws(x, average)
  if (summary) {
    d <- summarise_column(d, "adjustment", probs)
  }
  d
}

#' Posterior draws of the covariate adjustments of a BSCM
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param average \[`logical(1)`]\cr If `TRUE`, the adjustments are averaged
#'   over the treated units within each posterior draw.
#' @noRd
covariate_adjustment_draws <- function(x, average = FALSE) {
  X <- get_Xs(x)
  time <- get_time(x)
  predictors <- x$setup$beta_names
  T_total <- get_T_total(x)
  b <- rep(rvars_of(x, "beta"), each = T_total)
  dim(b) <- c(T_total, length(predictors))
  if (has_tv_coefs(x)) {
    tv_idx <- x$setup$tv_idx
    b[, tv_idx] <- b[, tv_idx] + rvars_of(x, "gamma")
  }
  d <- lapply(
    seq_len(get_N(x)),
    \(i) {
      dplyr::tibble(
        variable = rep(predictors, each = T_total),
        "{time}" := rep(get_times(x), times = length(predictors)),
        adjustment = c(b * covariate_delta(x, i, X))
      )
    }
  ) |>
    stats::setNames(get_treated(x)) |>
    dplyr::bind_rows(.id = get_unit(x))
  if (average && get_N(x) > 1L) {
    d <- d |>
      dplyr::summarise(
        adjustment = posterior::rvar_mean(.data$adjustment),
        .by = dplyr::all_of(c("variable", time))
      )
  }
  d
}
