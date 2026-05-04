#' Visualize donor weights
#'
#' `plot_weights()` visualizes posterior summaries of donor weights for a
#' fitted Bayesian synthetic control model. When applied to the output of
#' [leave_donor_out()], it instead visualizes how donor
#' weights or their ranks change across the leave-out runs.
#'
#' @param x \[`bscmfit` or `bscm_ldo`]\cr Output from [bscm()] or
#'   [leave_donor_out()].
#' @param point_estimate \[´character(1)`] Should the point estimate in weight
#'   plot correspond to posterior `"median"` (the default), or `"mean"`?
#' @param order \[`character()` or `NULL`]\cr Order of donors in y-axis. Either
#'   a vector of donor names, or `"descending"` / `"ascending"` to
#'   rank donors by posterior means of donor weights. Default is `"descending"`.
#'   Use `NULL` to use the original donor order. Ignored for `bscm_ldo` method
#'   which uses the ordering which was used in [leave_donor_out()].
#' @param coverage \[`numeric()`]\cr Coverages of posterior intervals of donor
#'   weights when `type = "weight"`. By default 50% and 95% posterior
#'   intervals are drawn.
#' @param type \[`character(1)`]\cr For `bscm_ldo` objects, plot posteriors
#'   of donor weights (`"weight"`) or donor weight ranks (`"rank"`).
#'   Default is `"weight"`, which is also only option for `bscmfit` objects.
#' @param linewidth \[`numeric(1)`]\cr Maximum line width used for intervals or
#'   donor rank trajectories. Default is `1`.
#' @param point_size \[`numeric(1)`]\cr Point size for point estimates or
#'   trajectory endpoints. Default is `2`.
#' @param reverse \[`logical(1)`]\cr Should the y-axis be reversed? Default is
#'   `NULL`, in which case y-axis is reversed except when
#'   `order  = "ascending"`, so in both ascending and descending ordering of
#'   donors leads to plot where donor with largest weight is on top.
#' @param ... Ignored.
#' @return A `ggplot` object, or a named list of `ggplot` objects for models
#'   with multiple treated units.
#' @export
plot_weights <- function(x, ...) {
  UseMethod("plot_weights", x)
}
#' @rdname plot_weights
#' @export
plot_weights.bscmfit <- function(
  x,
  point_estimate = "median",
  order = NULL,
  coverage = c(0.05, 0.95),
  linewidth = 1,
  point_size = 2,
  reverse = NULL,
  ...
) {
  stopifnot_(
    all(coverage > 0 & coverage < 1),
    "Argument {.arg coverage} must have values between 0 and 1."
  )
  point_estimate <- try_(match.arg(point_estimate, c("median", "mean")))
  stopifnot_(
    !inherits(point_estimate, "try-error"),
    "Argument {.arg point_estimate} must be either {.val median} or 
    {.val mean}."
  )
  stopifnot_(
    checkmate::test_number(linewidth, lower = 0, finite = TRUE),
    "Argument {.arg linewidth} must be a single positive number."
  )
  stopifnot_(
    checkmate::test_number(point_size, lower = 0, finite = TRUE),
    "Argument {.arg point_size} must be a single positive number."
  )
  stopifnot_(
    checkmate::test_flag(reverse, null.ok = TRUE),
    "Argument {.arg reverse} be single logical value or {.cls NULL}."
  )
  alpha <- (1 - sort(coverage)) / 2
  probs <- c(alpha, 1 - alpha)
  if (point_estimate == "median") {
    probs <- c(0.5, probs)
  }
  treated <- get_treated(x)
  N <- get_N(x)
  plots <- stats::setNames(vector("list", N), treated)
  unit <- get_unit(x)
  weights <- donor_weights(x, probs = probs) |>
    mutate(unit = .data[[unit]])
  donors <- order_donors(x, order, weights)
  if (is.null(reverse)) {
    reverse <- ifelse(order == "ascending", FALSE, TRUE)
  }
  for (i in treated) {
    plots[[i]] <- weights |>
      filter(.data$treated_unit == i) |>
      weight_plot(
        coverage,
        point_estimate,
        donors,
        linewidth,
        point_size,
        reverse = reverse
      )
  }
  if (length(plots) == 1L) plots[[1L]] else plots
}
#' @rdname plot_weights
#' @export
plot_weights.bscm_ldo <- function(
  x,
  type = "weight",
  point_estimate = "median",
  coverage = 0.95,
  linewidth = 1,
  point_size = 2,
  reverse = TRUE,
  ...
) {
  removed_donor <- treated_unit <- NULL
  type <- try_(match.arg(type, c("weight", "rank")))
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val weight} or {.val rank}."
  )
  point_estimate <- try_(match.arg(point_estimate, c("median", "mean")))
  stopifnot_(
    !inherits(point_estimate, "try-error"),
    "Argument {.arg point_estimate} must be either {.val median} or 
    {.val mean}."
  )
  stopifnot_(
    all(coverage > 0 & coverage < 1),
    "Argument {.arg coverage} must have values between 0 and 1."
  )
  stopifnot_(
    checkmate::test_number(linewidth, lower = 0, finite = TRUE),
    "Argument {.arg linewidth} must be a single positive number."
  )
  stopifnot_(
    checkmate::test_number(point_size, lower = 0, finite = TRUE),
    "Argument {.arg point_size} must be a single positive number."
  )
  stopifnot_(
    checkmate::test_flag(reverse, null.ok = TRUE),
    "Argument {.arg reverse} be single logical value or {.cls NULL}."
  )
  if (is.null(reverse)) {
    reverse <- ifelse(x$metadata$order == "ascending", FALSE, TRUE)
  }
  x$weights <- x$weights |> mutate(unit = .data[[x$metadata$unit]])
  donor_order <- x$metadata$donor_order
  treated <- x$metadata$treated
  plots <- stats::setNames(vector("list", length(treated)), treated)

  if (identical(type, "rank")) {
    keep <- c("removed_donor", "treated_unit", "unit", "mean")
    weights <- x$weights |>
      mutate(
        removed_donor = factor(removed_donor, levels = c("none", donor_order))
      ) |>
      select(all_of(keep)) |>
      filter(unit %in% donor_order) |>
      arrange(unit, removed_donor)

    for (i in treated) {
      plots[[i]] <- rank_plot(
        weights |> filter(.data$treated_unit == i),
        donor_order,
        linewidth,
        point_size,
        reverse
      )
    }
  } else {
    weights <- x$weights |> filter(unit %in% donor_order)
    for (i in treated) {
      plots[[i]] <- weight_plot_ldo(
        weights |> filter(treated_unit == i),
        coverage,
        point_estimate,
        donor_order,
        linewidth,
        point_size,
        reverse
      )
    }
  }
  if (length(plots) == 1L) plots[[1L]] else plots
}

rank_plot <- function(
  weights,
  donor_order,
  linewidth,
  point_size,
  reverse = TRUE
) {
  removed_donor <- x <- gap <- segment <- x_next <- rank_next <- NULL
  rank_data <- weights |>
    mutate(
      rank = rank(-mean, ties.method = "first"),
      .by = removed_donor
    ) |>
    mutate(
      x = as.integer(removed_donor),
      gap = (x - dplyr::lag(x, default = first(x))) != 1L,
      segment = cumsum(gap),
      .by = unit
    )

  breaks <- seq_along(donor_order)
  labels <- rank_data |>
    filter(removed_donor == "none") |>
    arrange(rank) |>
    pull(unit)
  rank_data <- rank_data |>
    mutate(unit = ordered(unit, levels = labels))

  point_data <- rank_data |>
    filter(
      row_number() == 1 | row_number() == n(),
      .by = c(unit, segment)
    )
  gap_data <- rank_data |>
    mutate(
      x_next = dplyr::lead(x),
      rank_next = dplyr::lead(rank),
      gap = (x_next - x) > 1,
      .by = unit
    ) |>
    filter(gap)

  if (reverse) {
    scale_y <- scale_y_reverse(
      NULL,
      labels = labels,
      breaks = breaks,
      expand = expansion(mult = 0.02, add = 0)
    )
  } else {
    scale_y <- scale_y_continuous(
      NULL,
      labels = labels,
      breaks = breaks,
      expand = expansion(mult = 0.02, add = 0)
    )
  }
  rank_data |>
    ggplot(aes(removed_donor, rank, colour = unit)) +
    geom_segment(
      data = gap_data,
      aes(
        xend = x_next,
        yend = rank_next
      ),
      linewidth = 0.5 * linewidth,
      alpha = 0.5
    ) +
    geom_line(
      aes(group = interaction(unit, segment)),
      linewidth = linewidth,
      alpha = 0.7
    ) +
    geom_point(
      data = point_data,
      size = point_size
    ) +
    scale_y +
    labs(x = "Removed donor") +
    scale_colour_viridis_d(guide = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}
weight_plot <- function(
  x,
  coverage,
  point_estimate,
  donors,
  linewidth,
  point_size,
  ldo_points = NULL,
  reverse = TRUE
) {
  point_col <- ifelse(point_estimate == "median", "q50", "mean")

  alpha <- (1 - coverage) / 2
  intervals <- sort(alpha[alpha != 0.5])
  lwr <- paste0("q", 100 * intervals)
  upr <- paste0("q", 100 * (1 - intervals))
  x <- x |> filter(!is.na(unit))
  n <- length(intervals)
  if (n == 1L) {
    alpha <- 1
  } else {
    z <- (seq_len(n) - 1) / (n - 1)
    linewidth <- (0.5 + 0.5 * z) * linewidth
    alpha <- 0.25 + 0.75 * z
  }
  if (reverse) {
    scale_y <- scale_y_reverse(
      NULL,
      expand = expansion(mult = 0.02, add = 0)
    )
  } else {
    scale_y <- scale_y_continuous(
      NULL,
      expand = expansion(mult = 0.02, add = 0)
    )
  }
  p <- ggplot(x, aes(.data[[point_col]], unit))
  if (!is.null(ldo_points)) {
    p <- p + ldo_points
  }
  for (i in seq_along(intervals)) {
    p <- p +
      geom_linerange(
        aes(xmin = .data[[lwr[i]]], xmax = .data[[upr[i]]]),
        linewidth = linewidth[i],
        alpha = alpha[i]
      )
  }
  p +
    geom_point(size = point_size) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.1),
      expand = expansion(mult = 0.02)
    ) +
    scale_y +
    labs(x = "Donor weight") +
    theme_bw()
}
weight_plot_ldo <- function(
  x,
  coverage,
  point_estimate,
  donors,
  linewidth,
  point_size,
  reverse
) {
  labels <- x |>
    filter(.data$removed_donor == "none") |>
    arrange(mean) |>
    pull(unit)
  x <- x |>
    mutate(unit = ordered(unit, levels = labels))

  base <- filter(x, .data$removed_donor == "none")
  ldo <- filter(x, .data$removed_donor != "none")
  check_weight_columns(x, c(if (point_estimate == "median") 0, coverage))

  ldo_points <- geom_point(
    data = ldo,
    aes(x = mean, y = unit),
    shape = 124,
    size = 1.5 * point_size,
    alpha = 0.7,
    colour = "grey50"
  )
  weight_plot(
    base,
    coverage,
    point_estimate,
    donors,
    linewidth,
    point_size,
    ldo_points,
    reverse
  )
}

check_weight_columns <- function(x, coverage) {
  alpha <- (1 - coverage) / 2
  expected_names <- paste0("q", 100 * c(alpha, 1 - alpha))
  missing <- setdiff(expected_names, names(x))
  stopifnot_(
    length(missing) == 0L,
    c(
      paste0(
        "Quantiles of posterior weights matching {.arg coverage} are not ",
        "present in the donor weights."
      ),
      i = "No columns with names {missing} in {.code weights}."
    )
  )
}
