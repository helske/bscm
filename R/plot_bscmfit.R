#' Visualize BSCM estimates
#'
#' A plot of the posterior mean and posterior interval of the treatment effect
#' and synthetic control over time for single treated unit models. For models
#' with multiple treated units, plots the average treatment effect over time
#' since treatment across all treated units.
#'
#' @param x \[`bscmfit`]\cr object.
#' @param unit \[`character()`]\cr Treated units to plot. If `NULL` (the
#'   default), all treated units are plotted.
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of
#' the posterior interval. Default is `c(0.025, 0.975)`.
#' @param ... Ignored
#' @aliases plot
#' @return A `ggplot` object when a single treated unit is plotted, and a
#'   named list of `ggplot` objects otherwise.
#' @export
#' @examples
#' plot(fit_single_treated, probs = c(0.05, 0.95))
plot.bscmfit <- function(x, unit = NULL, probs = c(0.025, 0.975), ...) {
  probs <- interval_probs(probs)
  outcome <- get_outcome(x)
  time <- get_time(x)
  unit_col <- get_unit(x)
  treatment <- get_treatment(x)
  treated <- get_treated(x)
  units <- check_units(x, unit)
  d_plot <- dplyr::bind_rows(
    "Treatment effect" = effects_draws(x) |>
      summarise_column("effect", probs, for_plots = TRUE),
    "Synthetic control" = sc_draws(x) |>
      summarise_column("y_rep", probs, for_plots = TRUE),
    .id = "type"
  )
  lookup <- stats::setNames(
    c(unit_col, time, treatment, outcome, paste0("q", 100 * probs)),
    c("unit", "time", "treatment", "mean", "ymin", "ymax")
  )
  d_plot <- d_plot |>
    dplyr::rename(dplyr::any_of(lookup))
  dy <- x$data |>
    dplyr::select(dplyr::all_of(c(unit_col, time, treatment, outcome))) |>
    dplyr::rename(dplyr::any_of(lookup)) |>
    dplyr::mutate(type = "Synthetic control", treatment = factor(treatment))

  start_t <- stats::setNames(get_times(x)[get_T_pre(x) + 1], treated)

  plots <- lapply(
    units,
    \(u) plot_bscm_unit(d_plot, dy, u, start_t[u])
  )
  if (length(plots) == 1L) {
    plots[[1L]]
  } else {
    stats::setNames(plots, units)
  }
}

#' Plot BSCM estimates for one treated unit
#' @noRd
plot_bscm_unit <- function(d, dy, id, start_t) {
  type <- time <- yintercept <- xintercept <- ymin <- ymax <- NULL
  dt <- data.frame(
    yintercept = 0,
    xintercept = start_t,
    type = "Treatment effect"
  )

  d |>
    dplyr::filter(unit == id) |>
    ggplot(aes(time, mean)) +
    geom_hline(
      data = dt,
      aes(yintercept = yintercept),
      linetype = "dashed",
      colour = "grey70"
    ) +
    geom_vline(
      data = dt,
      aes(xintercept = xintercept),
      linetype = "dashed",
      colour = "grey70"
    ) +
    geom_ribbon(
      aes(ymin = ymin, ymax = ymax, fill = type),
      alpha = 0.5
    ) +
    geom_point(
      data = dy |>
        dplyr::filter(unit == id),
      colour = "grey30"
    ) +
    geom_line(aes(colour = type)) +
    labs(x = time, y = NULL, subtitle = id) +
    scale_colour_manual(values = c("#DDAA33", "#0C7BDC")) +
    scale_fill_manual(values = c("#EECC66", "#77AADD")) +
    guides(fill = "none", colour = "none") +
    facet_grid(rows = vars(type), scales = "free_y") +
    theme_bw()
}
