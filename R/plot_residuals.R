#' Visualize BSCM residuals
#'
#' `plot_residuals()` plots the posterior mean and posterior interval of the
#' pre-treatment residuals over time (`type = "time"`), or the posterior
#' distribution of the residual autocorrelations up to a given number of lags
#' (`type = "autocorrelation"`), again only for pre-treatment time points.
#'
#' @param x \[`bscmfit`]\cr Object from [bscm()].
#' @param type \[`character(1)`]\cr Plot type. Either `"time"` (default) for
#'   a ribbon plot of residuals over time, or `"autocorrelation"` for posterior
#'   distributions of residual autocorrelations at each lag.
#' @param unit \[`character()`]\cr Treated units to plot. If `NULL` (the
#'   default), all treated units are plotted.
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of
#'   the posterior interval. Default is `c(0.025, 0.975)`.
#' @param max_lag \[`integer(1)`]\cr Maximum number of lags for
#'   `type = "autocorrelation"`. Default is `5L`.
#' @param ... Ignored.
#' @return A `ggplot` object when a single treated unit is plotted, and a
#'   named list of `ggplot` objects otherwise.
#' @seealso [residual_acf()] for the autocorrelations themselves, and
#'   [residuals.bscmfit()] for the residuals.
#' @export
plot_residuals <- function(x, ...) {
  UseMethod("plot_residuals", x)
}

#' @rdname plot_residuals
#' @export
#' @examples
#' plot_residuals(fit_single_treated)
#' plot_residuals(fit_single_treated, type = "autocorrelation", max_lag = 10L)
plot_residuals.bscmfit <- function(
  x,
  type = "time",
  unit = NULL,
  probs = c(0.025, 0.975),
  max_lag = 5L,
  ...
) {
  type <- try_(match.arg(type, c("time", "autocorrelation")))
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val time} or {.val autocorrelation}."
  )
  probs <- interval_probs(probs)
  units <- check_units(x, unit)
  time <- get_time(x)
  unit_col <- get_unit(x)
  qs <- paste0("q", 100 * probs)
  if (identical(type, "time")) {
    d <- residuals_draws(x) |>
      summarise_column("residuals", probs, for_plots = TRUE)
  } else {
    check_max_lag(max_lag)
    d <- residual_acf_draws(x, units, max_lag) |>
      summarise_column("ac", probs, for_plots = TRUE)
  }
  plots <- lapply(
    units,
    \(i) {
      d_i <- dplyr::filter(d, .data[[unit_col]] == .env$i)
      # the unit is named in the axis label only when the model has several
      id <- if (get_N(x) > 1L) i
      if (identical(type, "time")) {
        residual_time_plot(d_i, time, qs, id)
      } else {
        residual_acf_plot(d_i, qs, max_lag, id)
      }
    }
  )
  if (length(units) == 1L) {
    plots[[1L]]
  } else {
    stats::setNames(plots, units)
  }
}

#' Plot residuals of one treated unit over time
#' @noRd
residual_time_plot <- function(d, time, qs, id) {
  ggplot(d, aes(.data[[time]], .data$mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_ribbon(
      aes(ymin = .data[[qs[1]]], ymax = .data[[qs[2]]]),
      fill = "#77AADD",
      alpha = 0.25
    ) +
    geom_line(colour = "#0C7BDC") +
    labs(
      x = time,
      y = paste0("Residuals", if (!is.null(id)) paste0(" for ", id))
    ) +
    theme_bw()
}

#' Plot residual autocorrelations of one treated unit
#' @noRd
residual_acf_plot <- function(d, qs, max_lag, id) {
  ggplot(d, aes(.data$lag, .data$mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_ribbon(
      aes(ymin = .data[[qs[1]]], ymax = .data[[qs[2]]]),
      fill = "#77AADD",
      alpha = 0.25
    ) +
    geom_line(colour = "#0C7BDC") +
    geom_point(colour = "#0C7BDC") +
    scale_x_continuous("Lag", breaks = seq_len(max_lag)) +
    labs(
      y = paste0("Autocorrelation", if (!is.null(id)) paste0(" for ", id))
    ) +
    theme_bw()
}

#' @export
#' @rdname residual_acf
residual_acf <- function(x, ...) {
  UseMethod("residual_acf", x)
}
#' Autocorrelation function for residuals of BSCM
#'
#' Returns posterior draws (or summaries) of the autocorrelations of the
#' pre-treatment residuals up to lag `max_lag`, computed separately for each
#' treated unit.
#'
#' @inheritParams bscm_postprocessing
#' @param max_lag \[`integer(1)`]\cr Positive integer defining the maximum lag
#'   at which to compute the autocorrelations. The default is `5L`.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   posterior draws (`summary = FALSE`) in long format.
#' @rdname residual_acf
#' @aliases residual_acf
#' @seealso [plot_residuals()] for visualizing the autocorrelations.
#' @export
#' @examples
#' residual_acf(fit_single_treated, max_lag = 10)
residual_acf.bscmfit <- function(
  x,
  max_lag = 5L,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  check_flag(summary, "summary")
  check_max_lag(max_lag)
  probs <- sort_probs(probs)
  check_has_data(x, "x")
  d <- residual_acf_draws(x, get_treated(x), max_lag)
  if (summary) {
    d <- summarise_column(d, "ac", probs)
  }
  d
}

#' Posterior draws of the residual autocorrelations of a BSCM
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param units \[`character()`]\cr Treated units to compute for.
#' @param max_lag \[`integer(1)`]\cr Maximum lag.
#' @noRd
residual_acf_draws <- function(x, units, max_lag) {
  unit_col <- get_unit(x)
  d <- residuals_draws(x)
  lapply(
    units,
    \(i) {
      e <- d |>
        dplyr::filter(.data[[unit_col]] == .env$i) |>
        dplyr::pull("residuals") |>
        posterior::draws_of(with_chains = TRUE)
      n_iter <- dim(e)[1L]
      n_chains <- dim(e)[2L]
      acf_draws <- array(0, dim = c(n_iter, n_chains, max_lag))
      for (chain in seq_len(n_chains)) {
        for (s in seq_len(n_iter)) {
          ac <- stats::acf(
            e[s, chain, ],
            lag.max = max_lag,
            plot = FALSE
          )$acf
          # the first element is the lag zero autocorrelation
          acf_draws[s, chain, ] <- c(ac)[-1L]
        }
      }
      dplyr::tibble(
        lag = seq_len(max_lag),
        ac = posterior::rvar(acf_draws, with_chains = TRUE)
      )
    }
  ) |>
    stats::setNames(units) |>
    dplyr::bind_rows(.id = unit_col)
}
