#' Visualize BSCM treatment effects
#'
#' `plot_effects()` plots the posterior mean and posterior interval of the
#' treatment effect over time. For the output of
#' [leave_donor_out()] or [placebo_effects()], it additionally overlays the
#' posterior mean of the treatment effect from each leave-out or placebo fit.
#'
#' For models with multiple treated units, `plot_effects.bscmfit()` returns one
#' plot per treated unit. Supply a `unit` identifier to plot a single unit, or
#' `average = TRUE` to plot the average treatment effect over time since
#' treatment instead.
#'
#' @param x \[`bscmfit`, `bscm_ldo`, or `bscm_placebo_effects`]\cr Object from
#'   [bscm()], [leave_donor_out()], or [placebo_effects()].
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of
#'   the posterior interval. Default is `c(0.025, 0.975)`. For `bscm_ldo` and
#'   `bscm_placebo_effects` objects, defaults to the outermost probabilities
#'   used in [leave_donor_out()] or [placebo_effects()].
#' @param unit \[`character()`]\cr Treated units to plot. If `NULL` (the
#'   default), all treated units are plotted.
#' @param average \[`logical(1)`]\cr If `TRUE`, plots the average treatment
#'   effect over treated units at each time since treatment, instead of
#'   unit-specific effects. The default is `FALSE`.
#' @param ... Ignored.
#' @return A `ggplot` object when a single treated unit is plotted, and a
#'   named list of `ggplot` objects otherwise.
#' @export
plot_effects <- function(x, ...) {
  UseMethod("plot_effects", x)
}
#' @rdname plot_effects
#' @export
#' @examples
#' plot_effects(fit_single_treated)
plot_effects.bscmfit <- function(
  x,
  probs = c(0.025, 0.975),
  unit = NULL,
  average = FALSE,
  ...
) {
  probs <- interval_probs(probs)
  check_flag(average, "average")
  units <- check_units(x, unit)
  time <- get_time(x)
  unit_col <- get_unit(x)
  qs <- paste0("q", 100 * probs)
  if (average && get_N(x) > 1L) {
    d <- effects_draws(x, average = TRUE) |>
      summarise_column("effect", probs, for_plots = TRUE)
    return(average_effect_plot(d, time, qs))
  }
  d <- effects_draws(x) |>
    summarise_column("effect", probs, for_plots = TRUE)
  plots <- lapply(
    units,
    \(i) {
      d_i <- d |>
        dplyr::filter(.data[[unit_col]] == .env$i) |>
        dplyr::select(-dplyr::all_of(unit_col))
      unit_effect_plot(
        d_i,
        time,
        qs,
        start_t = d_i[[time]][get_T_pre(x)[i] + 1],
        id = if (get_N(x) > 1L) i
      )
    }
  )
  if (length(plots) == 1L) {
    plots[[1L]]
  } else {
    stats::setNames(plots, units)
  }
}

#' Plot the average treatment effect over time since treatment
#' @noRd
average_effect_plot <- function(d, time, qs) {
  ggplot(d, aes(.data[[time]], .data$mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_ribbon(
      aes(ymin = .data[[qs[1]]], ymax = .data[[qs[2]]]),
      fill = "#EECC66",
      alpha = 0.25
    ) +
    geom_line(colour = "#DDAA33") +
    labs(x = "Time since treatment", y = "Average treatment effect") +
    theme_bw()
}

#' Plot the treatment effect of one treated unit over time
#' @noRd
unit_effect_plot <- function(d, time, qs, start_t, id) {
  ggplot(d, aes(.data[[time]], .data$mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = start_t, linetype = "dashed", colour = "grey70") +
    geom_ribbon(
      aes(ymin = .data[[qs[1]]], ymax = .data[[qs[2]]]),
      fill = "#EECC66",
      alpha = 0.25
    ) +
    geom_line(colour = "#DDAA33") +
    labs(
      x = time,
      y = paste0("Treatment effect", if (!is.null(id)) paste0(" for ", id))
    ) +
    theme_bw()
}

#' @rdname plot_effects
#' @export
plot_effects.bscm_ldo <- function(x, probs = NULL, ...) {
  ymin <- ymax <- removed_donor <- NULL
  if (is.null(probs)) {
    probs <- range(x$metadata$probs)
  }
  probs <- interval_probs(probs)
  ymin_col <- paste0("q", 100 * min(probs))
  ymax_col <- paste0("q", 100 * max(probs))
  stopifnot_(
    ymin_col %in% names(x$effect) && ymax_col %in% names(x$effect),
    c(
      "Quantile columns matching {.arg probs} are not in the effect data.",
      i = paste0(
        "Re-run {.fun leave_donor_out} with {.arg probs} that include",
        "{.val {c(min(probs), max(probs))}}."
      )
    )
  )
  time <- x$metadata$time
  unit <- x$metadata$unit
  treated <- x$metadata$treated
  N <- length(treated)
  lookup <- stats::setNames(
    c(unit, time, ymin_col, ymax_col),
    c("unit", "time", "ymin", "ymax")
  )

  d_base <- x$effect |>
    dplyr::filter(removed_donor == "none") |>
    dplyr::rename(dplyr::any_of(lookup))

  d_ldo <- x$effect |>
    dplyr::filter(removed_donor != "none") |>
    dplyr::rename(dplyr::any_of(lookup))

  if (N > 1L) {
    ylab <- "Average treatment effect"
  } else {
    ylab <- "Treatment effect"
  }
  d_base |>
    ggplot(aes(time, mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_ribbon(
      aes(ymin = ymin, ymax = ymax),
      fill = "#EECC66",
      alpha = 0.25
    ) +
    geom_line(
      data = d_ldo,
      aes(group = removed_donor),
      colour = "grey50",
      alpha = 0.5
    ) +
    geom_line(colour = "#DDAA33") +
    labs(x = time, y = ylab) +
    theme_bw()
}
#' @rdname plot_effects
#' @export
plot_effects.bscm_placebo_effects <- function(x, probs = NULL, ...) {
  ymin <- ymax <- placebo <- NULL
  if (is.null(probs)) {
    probs <- range(x$metadata$probs)
  }
  probs <- interval_probs(probs)
  ymin_col <- paste0("q", 100 * min(probs))
  ymax_col <- paste0("q", 100 * max(probs))
  stopifnot_(
    ymin_col %in% names(x$effect) && ymax_col %in% names(x$effect),
    c(
      "Quantile columns matching {.arg probs} are not in the effect data.",
      i = paste0(
        "Re-run {.fun placebo_effects} with {.arg probs} that include ",
        "{.val {c(min(probs), max(probs))}}."
      )
    )
  )

  type <- x$metadata$type
  time <- x$metadata$setup$time
  unit <- x$metadata$setup$unit
  treated <- x$metadata$setup$treated

  lookup <- stats::setNames(
    c(time, ymin_col, ymax_col),
    c("time", "ymin", "ymax")
  )
  if (identical(type, "donor")) {
    d_base <- x$effect |>
      dplyr::filter(placebo == .env$treated) |>
      dplyr::rename(dplyr::any_of(lookup))
    d_placebo <- x$effect |>
      dplyr::filter(placebo != .env$treated) |>
      dplyr::rename(dplyr::any_of(lookup))
  } else {
    T_pre <- x$metadata$setup$T_pre
    d_base <- x$effect |>
      dplyr::filter(placebo == .env$T_pre) |>
      dplyr::rename(dplyr::any_of(lookup))
    d_placebo <- x$effect |>
      dplyr::filter(placebo != .env$T_pre) |>
      dplyr::rename(dplyr::any_of(lookup))
  }
  d_base |>
    ggplot(aes(time, mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_ribbon(
      aes(ymin = ymin, ymax = ymax),
      fill = "#EECC66",
      alpha = 0.25
    ) +
    geom_line(
      data = d_placebo,
      aes(group = placebo),
      colour = "grey50",
      alpha = 0.5
    ) +
    geom_line(colour = "#DDAA33") +
    labs(x = time, y = "Treatment effect") +
    theme_bw()
}
