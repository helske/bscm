#' Visualize BSCM residuals
#'
#' `plot_residuals()` plots the posterior mean and posterior interval of the
#' pre-treatment residuals over time (`type = "time"`), or the posterior
#' distribution of the residual autocorrelations up to a given number of lags
#' (`type = "autocorrelation"`), again only for pre-treatment time points.
#'
#' For models with multiple treated units, returns a named list of per-unit
#' plots.
#'
#' @param x \[`bscmfit`]\cr Object from [bscm()].
#' @param type \[`character(1)`]\cr Plot type. Either `"time"` (default) for
#'   a ribbon plot of residuals over time, or `"autocorrelation"` for posterior
#'   distributions of residual autocorrelations at each lag.
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of
#'   the posterior interval. Default is `c(0.025, 0.975)`.
#' @param max_lag \[`integer(1)`]\cr Maximum number of lags for
#'   `type = "autocorrelation"`. Default is `5L`.
#' @param ... Ignored.
#' @return A `ggplot` object, or a named list of `ggplot` objects when `x` is
#'   a `bscmfit` with multiple treated units.
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
  probs = c(0.025, 0.975),
  max_lag = 5L,
  ...
) {
  type <- try_(match.arg(type, c("time", "autocorrelation")))
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val time} or {.val autocorrelation}."
  )
  probs <- sort(probs)
  stopifnot_(
    length(probs) == 2L,
    "Argument {.arg probs} must be a {.cls numeric} vector of length 2 with 
    values between 0 and 1."
  )

  treated <- get_treated(x)
  time <- get_time(x)
  unit <- get_unit(x)
  N <- get_N(x)

  plots <- stats::setNames(vector("list", N), treated)

  lookup <- stats::setNames(
    c(time, paste0("q", 100 * probs)),
    c("time", "ymin", "ymax")
  )
  if (identical(type, "time")) {
    d <- residuals(x, probs = probs)
    d <- d |> dplyr::rename(dplyr::any_of(lookup))
    for (i in treated) {
      d_i <- if (N > 1L) dplyr::filter(d, .data[[unit]] == i) else d
      plots[[i]] <- d_i |>
        ggplot(aes(.data$time, .data$mean)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
        geom_ribbon(
          aes(ymin = .data$ymin, ymax = .data$ymax),
          fill = "#77AADD",
          alpha = 0.25
        ) +
        geom_line(colour = "#0C7BDC") +
        labs(
          x = time,
          y = paste0("Residuals", if (N > 1L) paste0(" for ", i))
        ) +
        theme_bw()
    }
  } else {
    d_raw <- residuals(x, summary = FALSE)
    for (i in treated) {
      d_i <- if (N > 1L) dplyr::filter(d_raw, .data[[unit]] == i) else d_raw
      d_ac <- residual_acf(d_i, max_lag, probs) |>
        dplyr::rename(dplyr::any_of(lookup))
      plots[[i]] <- d_ac |>
        ggplot(aes(.data$lag, .data$mean)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
        geom_ribbon(
          aes(ymin = .data$ymin, ymax = .data$ymax),
          fill = "#77AADD",
          alpha = 0.25
        ) +
        geom_line(colour = "#0C7BDC") +
        geom_point(colour = "#0C7BDC") +
        scale_x_continuous("Lag", breaks = seq_len(max_lag)) +
        labs(
          y = paste0("Autocorrelation", if (N > 1L) paste0(" for ", i))
        ) +
        theme_bw()
    }
  }
  if (N == 1L) plots[[1L]] else plots
}

#' Autocorrelation function for residuals of BSCM
#'
#' @param x Posterior draws from [residuals.bscmfit()] with `summary = FALSE`,
#'   or `bscmfit` object.
#' @param max_lag  \[`integer(1)`]\cr Positive integer defining the maximum lag
#' at which to compute the autocorrelations. The default is `5L`.
#' @inheritParams plot_residuals
#' @return A data frame.
#' @export
#' @examples
#' residual_acf(fit_single_treated, max_lag = 10)
residual_acf <- function(x, max_lag = 5L, probs = c(0.025, 0.975)) {
  stopifnot_(
    checkmate::test_count(max_lag),
    "Argument {.arg max_lag} must be a single positive integer."
  )
  probs <- sort(probs)
  if (inherits(x, "bscmfit")) {
    x <- residuals(x, summary = FALSE)
  }
  e <- t(posterior::as_draws_matrix(x$residuals))
  n_draws <- ncol(e)
  acf_mat <- matrix(0, max_lag, n_draws)
  for (s in seq_len(n_draws)) {
    ac <- c(stats::acf(e[, s], lag.max = max_lag, plot = FALSE)$acf)
    acf_mat[, s] <- ac[-1L]
  }
  d <- dplyr::tibble(
    lag = seq_len(max_lag),
    ac = posterior::rvar(t(acf_mat))
  )
  d |>
    dplyr::mutate(summarise_with_probs(.data$ac, probs, for_plots = TRUE)) |>
    dplyr::select(-"ac", -"variable")
}
