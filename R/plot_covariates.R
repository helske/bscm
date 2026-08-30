#' @export
#' @rdname plot_covariate_imbalance
plot_covariate_imbalance <- function(x, ...) {
  UseMethod("plot_covariate_imbalance", x)
}
#' Plot the covariate imbalance of a Bayesian synthetic control model
#'
#' Plots the posterior mean and interval of the covariate imbalances over
#' time. See [covariate_imbalance()] for the definition of the imbalance.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param unit \[`character()`]\cr Treated units to plot. If `NULL` (the
#'   default), all treated units are plotted.
#' @param average \[`logical(1)`]\cr If `TRUE`, plots the imbalance averaged
#'   over the treated units instead of unit-specific imbalances. The default
#'   is `FALSE`.
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of
#'   the posterior interval. Default is `c(0.025, 0.975)`.
#' @param alpha \[`numeric(1)`]\cr Opacity of the posterior interval ribbon.
#'   Default is `0.25`.
#' @param ... Ignored.
#' @return A `ggplot` object when a single treated unit is plotted, and a
#'   named list of `ggplot` objects otherwise.
#' @rdname plot_covariate_imbalance
#' @aliases plot_covariate_imbalance
#' @seealso [covariate_imbalance()].
#' @export
#' @examples
#' plot_covariate_imbalance(fit_single_treated, probs = c(0.05, 0.95))
plot_covariate_imbalance.bscmfit <- function(
  x,
  unit = NULL,
  average = FALSE,
  probs = c(0.025, 0.975),
  alpha = 0.25,
  ...
) {
  probs <- interval_probs(probs)
  check_flag(average, "average")
  check_alpha(alpha)
  check_has_predictors(x)
  units <- check_units(x, unit)
  time <- get_time(x)
  qs <- paste0("q", 100 * probs)
  if (average && get_N(x) > 1L) {
    d <- covariate_imbalance_draws(x, average = TRUE) |>
      summarise_column("imbalance", probs, for_plots = TRUE)
    return(imbalance_plot(d, time, qs, alpha, id = NULL))
  }
  d <- covariate_imbalance_draws(x) |>
    summarise_column("imbalance", probs, for_plots = TRUE)
  unit_col <- get_unit(x)
  plots <- lapply(
    units,
    \(i) {
      d_i <- dplyr::filter(d, .data[[unit_col]] == .env$i)
      imbalance_plot(d_i, time, qs, alpha, id = if (get_N(x) > 1L) i)
    }
  )
  if (length(plots) == 1L) {
    plots[[1L]]
  } else {
    stats::setNames(plots, units)
  }
}

#' Plot the covariate imbalance of one treated unit
#' @noRd
imbalance_plot <- function(d, time, qs, alpha, id) {
  ggplot(d, aes(.data[[time]], .data$mean)) +
    geom_ribbon(
      aes(ymin = .data[[qs[1]]], ymax = .data[[qs[2]]]),
      alpha = alpha
    ) +
    geom_line() +
    labs(
      x = time,
      y = paste0("Covariate imbalance", if (!is.null(id)) paste0(" for ", id))
    ) +
    theme_bw()
}

#' @export
#' @rdname plot_covariate_adjustment
plot_covariate_adjustment <- function(x, ...) {
  UseMethod("plot_covariate_adjustment", x)
}
#' Plot the covariate adjustments of a Bayesian synthetic control model
#'
#' Plots the posterior mean and interval of the covariate adjustments over
#' time, with one panel per predictor. See [covariate_adjustment()] for the
#' definition of the adjustment.
#'
#' @inheritParams plot_covariate_imbalance
#' @param average \[`logical(1)`]\cr If `TRUE`, plots the adjustments averaged
#'   over the treated units instead of unit-specific adjustments. The default
#'   is `FALSE`.
#' @param alpha \[`numeric(1)`]\cr Opacity of the posterior interval ribbon.
#'   Default is `0.5`.
#' @param scales \[`character(1)`]\cr Should the panel scales be free
#'   (`"free_y"`, the default) or fixed? Passed to [ggplot2::facet_wrap()].
#' @param ... Optional arguments passed to [ggplot2::facet_wrap()].
#' @return A `ggplot` object when a single treated unit is plotted, and a
#'   named list of `ggplot` objects otherwise.
#' @rdname plot_covariate_adjustment
#' @aliases plot_covariate_adjustment
#' @seealso [covariate_adjustment()].
#' @export
#' @examples
#' plot_covariate_adjustment(fit_single_treated)
plot_covariate_adjustment.bscmfit <- function(
  x,
  unit = NULL,
  average = FALSE,
  probs = c(0.025, 0.975),
  alpha = 0.5,
  scales = "free_y",
  ...
) {
  probs <- interval_probs(probs)
  check_flag(average, "average")
  check_alpha(alpha)
  check_has_predictors(x)
  units <- check_units(x, unit)
  time <- get_time(x)
  qs <- paste0("q", 100 * probs)
  if (average && get_N(x) > 1L) {
    d <- covariate_adjustment_draws(x, average = TRUE) |>
      summarise_column("adjustment", probs, for_plots = TRUE)
    return(adjustment_plot(d, time, qs, alpha, scales, id = NULL, ...))
  }
  d <- covariate_adjustment_draws(x) |>
    summarise_column("adjustment", probs, for_plots = TRUE)
  unit_col <- get_unit(x)
  plots <- lapply(
    units,
    \(i) {
      d_i <- dplyr::filter(d, .data[[unit_col]] == .env$i)
      adjustment_plot(
        d_i,
        time,
        qs,
        alpha,
        scales,
        id = if (get_N(x) > 1L) i,
        ...
      )
    }
  )
  if (length(plots) == 1L) {
    plots[[1L]]
  } else {
    stats::setNames(plots, units)
  }
}

#' Plot the covariate adjustments of one treated unit
#' @noRd
adjustment_plot <- function(d, time, qs, alpha, scales, id, ...) {
  ggplot(d, aes(.data[[time]], .data$mean)) +
    geom_ribbon(
      aes(ymin = .data[[qs[1]]], ymax = .data[[qs[2]]]),
      fill = "#77AADD",
      alpha = alpha
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(colour = "#0C7BDC") +
    scale_x_continuous(limits = range(d[[time]])) +
    facet_wrap(~variable, scales = scales, ...) +
    labs(
      x = time,
      y = paste0("Covariate adjustment", if (!is.null(id)) paste0(" for ", id))
    ) +
    theme_bw()
}
