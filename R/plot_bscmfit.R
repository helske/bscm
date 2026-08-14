#' Visualize BSCM estimates
#'
#' A plot of the posterior mean and posterior interval of the treatment effect
#' and synthetic control over time for single treated unit models. For models
#' with multiple treated units, plots the average treatment effect over time
#' since treatment across all treated units.
#'
#' @param x \[`bscmfit`]\cr object.
#' @param unit \[`character(1)`]\cr For models with
#'   multiple treated units, the name of the treated unit to 
#'   plot. If `NULL` (the default), all treated units are plotted sequentially.
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of
#' the posterior interval. Default is `c(0.025, 0.975)`.
#' @param ... Ignored
#' @aliases plot
#' @return A `ggplot` object or a list of `ggplot` objects in case of multiple 
#'   treated units when `unit` is `NULL`.
#' @export
#' @examples
#' plot(fit_single_treated, probs = c(0.05, 0.95))
plot.bscmfit <- function(x, unit = NULL, probs = c(0.025, 0.975), ...) {
  
  stopifnot_(
    checkmate::test_numeric(
      probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      len = 2L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector of length 2 with 
    values between 0 and 1."
  )
  outcome <- get_outcome(x) 
  time <- get_time(x) 
  unit_col <- get_unit(x) 
  treatment <- get_treatment(x)
  treated <- get_treated(x)
  if (!is.null(unit)) {
    stopifnot_(
      unit %in% treated,
      "Argument {.arg unit} must be one of the treated units:
      {.val {treated}}."
    )
  }
  d_plot <- dplyr::bind_rows(
    "Treatment effect" = treatment_effect(
      x, average = FALSE, probs = probs, for_plots = TRUE
    ),
    "Synthetic control" = synthetic_control(
      x, probs = probs, for_plots = TRUE
    ),
    .id = "type"
  ) 
  N <- get_N(x)
  if (N == 1L) unit <- treated[1]
  if (is.null(d_plot[[unit_col]])) {
    d_plot[[unit_col]] <- treated[1]
  }
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
 
  if (!is.null(unit)) {
    return(plot_bscm_unit(d_plot, dy, unit, start_t[unit]))
  }
  p <- lapply(
    treated,
    \(u) plot_bscm_unit(d_plot, dy, u, start_t[u])
  )
  names(p) <- treated
  p
}

#' Plot BSCM estimates for one treated unit
#' @noRd 
plot_bscm_unit <- function(d, dy, id, start_t) { 
  type <- time <- yintercept <- xintercept <- ymin <- ymax <- NULL
  dt <- data.frame(
    yintercept = 0,
    xintercept = start_t,
    type = "Treatment effect")
  
  d |> 
    dplyr::filter(unit == id) |> 
    ggplot(aes(time, mean)) +
    geom_hline(
      data = dt, aes(yintercept = yintercept),
      linetype = "dashed", colour = "grey70"
    ) +
    geom_vline(
      data = dt, aes(xintercept = xintercept),
      linetype = "dashed", colour = "grey70"
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
