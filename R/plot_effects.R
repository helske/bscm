#' Visualize BSCM treatment effects
#'
#' `plot_effects()` plots the posterior mean and posterior interval of the
#' treatment effect over time. For the output of
#' [leave_donor_out()] or [placebo_effects()], it additionally overlays the
#' posterior mean of the treatment effect from each leave-out or placebo fit.
#'
#' For models with multiple treated units, `plot_effects.bscmfit()` defaults
#' to plotting the average treatment effect over time since treatment. Supply a
#' `unit` identifier to plot the effect for a specific treated unit instead.
#'
#' @param x \[`bscmfit`, `bscm_ldo`, or `bscm_placebo_effects`]\cr Object from
#'   [bscm()], [leave_donor_out()], or [placebo_effects()].
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of
#'   the posterior interval. Default is `c(0.025, 0.975)`. For `bscm_ldo` and
#'   `bscm_placebo_effects` objects, defaults to the outermost probabilities
#'   used in [leave_donor_out()] or [placebo_effects()].
#' @param unit \[`character(1)` or `numeric(1)` or `NULL`]\cr For models with
#'   multiple treated units, the identifier of the specific treated unit to
#'   plot. If `NULL` (the default), the average treatment effect across all
#'   treated units is plotted.
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @export
plot_effects <- function(x, ...) {
  UseMethod("plot_effects", x)
}
#' @rdname plot_effects
#' @export
#' @examples
#' plot_effects(fit_single_treated)
plot_effects.bscmfit <- function(x, probs = c(0.025, 0.975), unit = NULL, ...) {
  ymin <- ymax <- NULL
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
  treated <- get_treated(x)
  time <- get_time(x)
  unit_col <- get_unit(x)
  N <- get_N(x)
  qcols <- paste0("q", 100 * sort(probs))

  if (!is.null(unit)) {
    stopifnot_(
      unit %in% treated,
      "Argument {.arg unit} must be one of the treated units:
      {.val {treated}}."
    )
  }

  if (is.null(unit) && N > 1L) {
    d <- treatment_effect(x, average = TRUE, probs = probs, for_plots = TRUE)
    d |>
      dplyr::rename(
        ymin = dplyr::all_of(qcols[[1L]]),
        ymax = dplyr::all_of(qcols[[2L]])
      ) |>
      ggplot(aes(.data[[time]], mean)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
      geom_ribbon(
        aes(ymin = ymin, ymax = ymax),
        fill = "#EECC66",
        alpha = 0.25
      ) +
      geom_line(colour = "#DDAA33") +
      labs(x = "Time since treatment", y = "Average treatment effect") +
      theme_bw()
  } else {
    d <- treatment_effect(x, average = FALSE, probs = probs, for_plots = TRUE)
    if (N > 1L) {
      d <- d |>
        dplyr::filter(.data[[unit_col]] == .env$unit) |>
        dplyr::select(-dplyr::all_of(unit_col))
    } else {
      unit <- treated[1]
    }
    start_t <- d[[time]][get_T_pre(x)[unit] + 1]
    d |>
      dplyr::rename(
        ymin = dplyr::all_of(qcols[[1L]]),
        ymax = dplyr::all_of(qcols[[2L]])
      ) |>
      ggplot(aes(.data[[time]], mean)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
      geom_vline(xintercept = start_t, linetype = "dashed", colour = "grey70") +
      geom_ribbon(
        aes(ymin = ymin, ymax = ymax),
        fill = "#EECC66",
        alpha = 0.25
      ) +
      geom_line(colour = "#DDAA33") +
      labs(
        x = time,
        y = paste0("Treatment effect", if (N > 1L) paste0(" for ", unit))
      ) +
      theme_bw()
  }
}
#' @rdname plot_effects
#' @export
plot_effects.bscm_ldo <- function(x, probs = NULL, ...) {
  ymin <- ymax <- removed_donor <- NULL
  if (is.null(probs)) {
    probs <- range(x$metadata$probs)
  }
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
