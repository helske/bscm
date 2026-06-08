#' Visualize BSCM estimates
#'
#' A plot of the posterior mean and posterior interval of the treatment effect
#' and synthetic control over time.
#'
#' @param x \[`bscmfit`]\cr object.
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of
#' the posterior interval. Default is `c(0.025, 0.975)`.
#' @param ... Ignored
#' @aliases plot
#' @return A `ggplot` object
#' @export
#' @examples
#' plot(fit_single_treated, probs = c(0.05, 0.95))
plot.bscmfit <- function(x, probs = c(0.025, 0.975), ...) {
  type <- yintercept <- ymin <- ymax <- NULL
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
  treated <- get_treated(x)
  time <- get_time(x)
  unit <- get_unit(x)
  treatment <- get_treatment(x)

  d_effects <- treatment_effect(
    x,
    type = "time",
    average = FALSE,
    probs = probs,
    for_plots = TRUE
  )
  d_synth <- synthetic_control(x, probs = probs, for_plots = TRUE)
  d_effects$type <- "Treatment effect"
  d_synth$type <- "Synthetic control"
  lookup <- stats::setNames(
    c(unit, time, treatment, outcome, paste0("q", 100 * probs)),
    c("unit", "time", "treatment", "mean", "ymin", "ymax")
  )
  d_plot <- dplyr::bind_rows(d_effects, d_synth) |>
    dplyr::rename(dplyr::any_of(lookup))
  dy <- x$data |>
    dplyr::filter(.data[[unit]] %in% .env$treated) |>
    dplyr::select(dplyr::all_of(c(unit, time, treatment, outcome))) |>
    dplyr::mutate(type = "Synthetic control", treatment = factor(treatment)) |>
    dplyr::rename(dplyr::any_of(lookup))

  wrap <- facet_grid(
    rows = vars(type),
    cols = if (get_N(x) > 1) vars(unit),
    scales = "free_y"
  )
  d_plot |>
    ggplot(aes(time, mean)) +
    geom_ribbon(
      aes(ymin = ymin, ymax = ymax, fill = type),
      alpha = 0.5
    ) +
    geom_line(aes(colour = type)) +
    geom_point(data = dy, aes(shape = treatment), colour = "grey30") +
    labs(x = time, y = NULL) +
    scale_colour_manual(values = c("#DDAA33", "#0C7BDC")) +
    scale_fill_manual(values = c("#EECC66", "#77AADD")) +
    guides(fill = "none", colour = "none", shape = "none") +
    wrap +
    theme_bw()
}
