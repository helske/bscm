#' Visualize regression coefficients of Bayesian synthetic control model
#'
#' Plots posterior means and posterior intervals for the model coefficients.
#' For time-constant coefficients (`type = "fixed"`), a horizontal point-range
#' plot is produced. For time-varying coefficients (`type = "varying"`), a
#' ribbon-and-line plot over time is produced, faceted by parameter.
#'
#' @param x \[`bscmfit`]\cr Output from [bscm()].
#' @param type \[`character(1)`]\cr Type of coefficients to plot. Either
#'   `"fixed"` for time-constant coefficients or `"varying"` for
#'   time-varying coefficients. Default is `"varying"` if the model contains
#'   varying coefficients.
#' @param combine \[`logical(1)`]\cr If `TRUE` and `type = "varying"`, plot
#'   the total effect (beta + gamma) instead of gamma alone.
#'   Ignored when `type = "fixed"`. Default is `TRUE`.
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
  combine = TRUE,
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
  check_flag(combine, "combine")
  if (type == "fixed") {
    plot_coefs_fixed(x, probs = probs)
  } else {
    stopifnot_(
      has_tv_coefs(x),
      "The model does not contain any time-varying predictors."
    )
    plot_coefs_varying(
      x,
      probs = probs,
      combine = combine,
      alpha = alpha,
      scales = scales
    )
  }
}

#' Plot Time-constant Coefficients of a \pkg{bscm} Model
#' @noRd
plot_coefs_fixed <- function(x, probs) {
  d <- coef_draws(x, "beta") |>
    summarise_column("beta", probs, for_plots = TRUE) |>
    dplyr::rename(variable = "parameter")
  qs <- paste0("q", 100 * probs)
  d$variable <- factor(d$variable, levels = d$variable)
  ggplot(d, aes(.data$mean, .data$variable)) +
    geom_linerange(
      aes(xmin = .data[[qs[1]]], xmax = .data[[qs[2]]])
    ) +
    geom_point() +
    labs(x = "Value", y = "Parameter") +
    theme_bw()
}


#' Plot Time-varying Coefficients of a \pkg{bscm} Model
#'
#' @noRd
plot_coefs_varying <- function(x, probs, combine, alpha, scales) {
  time_var <- get_time(x)
  qs <- paste0("q", 100 * probs)

  if (!combine) {
    d <- coef_draws(x, "gamma") |>
      summarise_column("gamma", probs, for_plots = TRUE) |>
      dplyr::rename(variable = "parameter")
  } else {
    d_beta <- coef_draws(x, "beta") |>
      dplyr::mutate(parameter = sub("^beta_", "", .data$parameter)) |>
      dplyr::filter(.data$parameter %in% x$setup$gamma_names)
    d_gamma <- coef_draws(x, "gamma") |>
      dplyr::mutate(parameter = sub("^gamma_", "", .data$parameter))
    d <- d_gamma |>
      dplyr::left_join(d_beta, by = "parameter") |>
      dplyr::mutate(
        total = .data$gamma + .data$beta,
        summarise_with_probs(.data$total, probs, for_plots = TRUE)
      ) |>
      dplyr::select(
        variable = "parameter",
        "mean",
        dplyr::all_of(c(time_var, qs))
      )
  }
  d |>
    ggplot(aes(.data[[time_var]], mean)) +
    geom_ribbon(
      aes(ymin = .data[[qs[1]]], ymax = .data[[qs[2]]]),
      fill = "#77AADD",
      alpha = alpha
    ) +
    geom_line(colour = "#0C7BDC") +
    scale_x_continuous(limits = range(d[[time_var]])) +
    facet_wrap(~variable, scales = scales) +
    labs(x = time_var, y = "Value") +
    theme_bw()
}
