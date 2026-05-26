#' @export
#' @rdname treatment_effect
treatment_effect <- function(x, ...) {
  UseMethod("treatment_effect", x)
}
#' Treatment effect estimates of a Bayesian synthetic control model
#'
#' @inheritParams rmse.bscmfit
#' @param type \[`character(1)`]\cr Type of treatment effect to compute.
#'   `"time"` (the default) returns effects at each time point.
#'   `"average"` returns the average pre- and post-treatment effects.
#'   `"cumulative"` returns the cumulative average effects over the
#'   post-treatment period.
#' @param average \[`logical(1)`]\cr If `TRUE` (the default), returns the
#' average effects over treated units in case of
#' multiple treated units. If `FALSE`, unit-specific effects are returned.
#' Currently not applicable to `type = "cumulative"`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or
#'   posterior samples (`summary = FALSE`) in long format.
#' @rdname treatment_effect
#' @aliases treatment_effect
#' @export
#' @examples
#' fit <- bscm(y ~ 1, data = single_treated, treatment = "treatment",
#'  chains = 1, control = list(adapt_delta = 0.8), refresh = 0)
#' treatment_effect(fit) |> tail()
treatment_effect.bscmfit <- function(
  x,
  type = "time",
  average = TRUE,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  probs <- sort_probs(probs)
  stopifnot_(
    checkmate::test_flag(summary),
    "Argument {.arg summary} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_flag(average),
    "Argument {.arg average} must be a single {.cls logical} value."
  )
  type <- try_(match.arg(type, c("time", "average", "cumulative")))
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val time}, {.val average} or 
    {.val cumulative}."
  )
  stopifnot_(
    !is.null(x$data),
    "The model fit {.arg x} does not contain the original data. You probably
    used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  time <- get_time(x)
  times <- get_times(x)
  T_pre <- get_T_pre(x)
  T_total <- get_T_total(x)
  N <- get_N(x)
  unit <- get_unit(x)
  treated <- get_treated(x)

  for_plots <- list(...)$for_plots %||% FALSE
  y_rep <- as_draws_rvars(as_draws(x, "y_rep"))$y_rep
  effect <- get_stan_y(x) - y_rep

  if (type == "average") {
    pre <- vector("list", N)
    post <- vector("list", N)
    for (i in seq_len(N)) {
      T_ <- T_pre[treated[i]]
      pre[[i]] <- rvar_mean(effect[seq_len(T_), i])
      post[[i]] <- rvar_mean(effect[seq.int(T_ + 1L, T_total), i])
    }
    if (average && N > 1) {
      vars <- c(
        "Average pre-treatment effect",
        "Average post-treatment effect"
      )
      values <- lapply(
        list(do.call(c, pre), do.call(c, post)),
        rvar_mean
      )
      out <- format_posterior_output(
        values,
        summary = summary,
        probs = probs,
        variable = vars,
        for_plots = for_plots
      )
    } else {
      values <- c(pre, post)
      vars <- rep(c("Pre-treatment effect", "Post-treatment effect"), each = N)
      units <- rep(treated, times = 2)
      out <- format_posterior_output(
        values,
        summary = summary,
        probs = probs,
        variable = vars,
        for_plots = for_plots
      ) |>
        add_output_column(
          name = unit,
          values = units,
          summary = summary
        )
    }
  }
  if (type == "cumulative") {
    cumavg <- vector("list", N)
    for (i in seq_len(N)) {
      T_ <- T_pre[treated[i]]
      cumavg[[i]] <- rvar_apply(
        effect[seq.int(T_ + 1L, T_total), i],
        2,
        cumsum
      )
    }
    out <- lapply(seq_len(N), \(i) {
      T_ <- T_pre[treated[i]]
      times_i <- times[seq.int(T_ + 1L, T_total)]
      format_posterior_output(
        cumavg[[i]],
        summary = summary,
        probs = probs,
        variable = "Cumulative post-treatment effect",
        for_plots = for_plots
      ) |>
        add_output_column(
          name = time,
          values = times_i,
          summary = summary
        ) |>
        add_output_column(
          name = unit,
          values = treated[i],
          summary = summary
        )
    }) |>
      bind_rows()
  }
  if (type == "time") {
    if (average && N > 1) {
      value <- effect |>
        rvar_apply(1, rvar_mean)
      out <- format_posterior_output(
        value,
        summary = summary,
        probs = probs,
        variable = "Average treatment effect",
        for_plots = for_plots
      ) |>
        add_output_column(
          name = time,
          values = times,
          summary = summary
        )
    } else {
      units <- rep(treated, each = length(times))
      times2 <- rep(times, times = length(treated))
      out <- format_posterior_output(
        effect,
        summary = summary,
        probs = probs,
        variable = "Treatment effect",
        for_plots = for_plots
      ) |>
        add_output_column(
          name = time,
          values = times2,
          summary = summary
        ) |>
        add_output_column(
          name = unit,
          values = units,
          summary = summary
        )
    }
  }
  out
}
