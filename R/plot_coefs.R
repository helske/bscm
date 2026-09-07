#' Visualize regression coefficients of Bayesian synthetic control model
#'
#' Plots posterior means and posterior intervals for the coefficients of the
#' predictors. For time-constant coefficients (`type = "fixed"`), a horizontal
#' point-range plot is produced. For time-varying coefficients
#' (`type = "varying"`), a ribbon-and-line plot of the total effect
#' `beta + gamma_t` over time is produced, faceted by parameter.
#'
#' @param x \[`bscmfit`]\cr Output from [bscm()].
#' @param type \[`character(1)`]\cr Type of coefficients to plot. Either
#'   `"fixed"` for time-constant coefficients or `"varying"` for
#'   time-varying coefficients. Default is `"varying"` if the model contains
#'   varying coefficients.
#' @param alpha \[`numeric(1)`]\cr Opacity of the credible-interval ribbon
#'   (used only when `type = "varying"`). Default is `0.5`.
#' @param scales \[`character(1)`]\cr Passed to [ggplot2::facet_wrap()]
#'   (used only when `type = "varying"`). Default is `"free_y"`.
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of
#'   the posterior interval. Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @export
plot_coefs <- function(x, ...) {
  UseMethod("plot_coefs", x)
}

#' @rdname plot_coefs
#' @export
#' @examples
#' plot_coefs(fit_single_treated)
plot_coefs.bscmfit <- function(
  x,
  type = NULL,
  probs = c(0.025, 0.975),
  alpha = 0.5,
  scales = "free_y",
  ...
) {
  stopifnot_(
    has_predictors(x),
    "The model does not contain any predictors."
  )

  if (is.null(type)) {
    type <- if (has_tv_coefs(x)) "varying" else "fixed"
  } else {
    type <- try_(match.arg(type, c("fixed", "varying")))
    stopifnot_(
      !inherits(type, "try-error"),
      "Argument {.arg type} must be either {.val fixed} or {.val varying}."
    )
  }

  probs <- interval_probs(probs)
  if (type == "fixed") {
    stopifnot_(
      length(fixed_idx(x)) > 0L,
      "The model does not contain any time-constant coefficients."
    )
    plot_coefs_fixed(x, probs = probs)
  } else {
    stopifnot_(
      has_tv_coefs(x),
      "The model does not contain any time-varying predictors."
    )
    plot_coefs_varying(
      x,
      probs = probs,
      alpha = alpha,
      scales = scales
    )
  }
}

#' Plot Time-constant Coefficients of a \pkg{bscm} Model
#' @noRd
plot_coefs_fixed <- function(x, probs) {
  d <- coef_draws(x, "fixed") |>
    summarise_column("coef", probs, diagnostics = FALSE) |>
    dplyr::mutate(
      parameter = factor(.data$parameter, levels = unique(.data$parameter))
    )
  qs <- paste0("q", 100 * probs)
  ggplot(d, aes(.data$mean, .data$parameter)) +
    geom_linerange(
      aes(xmin = .data[[qs[1]]], xmax = .data[[qs[2]]])
    ) +
    geom_point() +
    labs(x = "Value", y = "Parameter") +
    theme_bw()
}


#' Plot Time-varying Coefficients of a \pkg{bscm} Model
#'
#' Plots the total effect of each time-varying predictor, `beta + gamma_t`.
#'
#' @noRd
plot_coefs_varying <- function(x, probs, alpha, scales) {
  time_var <- get_time(x)
  qs <- paste0("q", 100 * probs)
  d <- coef_draws(x, "varying") |>
    summarise_column("coef", probs, diagnostics = FALSE) |> 
    dplyr::mutate(
      parameter = factor(.data$parameter, levels = unique(.data$parameter))
    )
  ggplot(d, aes(.data[[time_var]], .data$mean)) +
    geom_ribbon(
      aes(ymin = .data[[qs[1]]], ymax = .data[[qs[2]]]),
      fill = "#77AADD",
      alpha = alpha
    ) +
    geom_line(colour = "#0C7BDC") +
    scale_x_continuous(limits = range(d[[time_var]])) +
    facet_wrap(~ parameter, scales = scales) +
    labs(x = time_var, y = "Value") +
    theme_bw()
}
