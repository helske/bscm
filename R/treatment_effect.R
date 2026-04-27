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
#' @return A `data.frame` of posterior summaries of the treatment effects.
#' @rdname treatment_effect
#' @aliases treatment_effect
#' @export
treatment_effect.bscmfit <- function(x,
                                     type = "time",
                                     average = TRUE, 
                                     probs = c(0.025, 0.975), ...) {
  
  probs <- sort_probs(probs)
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
      out <- lapply(
        list(do.call(c, pre), do.call(c, post)), rvar_mean
      ) |> 
        summarise_with_probs(probs, for_plots) |>
        mutate(
          variable = c(
            "Average pre-treatment effect", "Average post-treatment effect"
          )
        )
    } else {
      out <- c(pre, post) |> 
        summarise_with_probs(probs, for_plots) |>
        mutate(
          "{unit}" := rep(treated, times = 2), .before = 1L
        ) |>
        mutate(
          variable = rep(
            c("Pre-treatment effect", "Post-treatment effect"), each = N
          )
        )
    }
  }
  if (type == "cumulative") {
    cumavg <-  vector("list", N)
    for (i in seq_len(N)) {
      T_ <- T_pre[treated[i]]
      cumavg[[i]] <- rvar_apply(
        effect[seq.int(T_ + 1L, T_total), i], 2, cumsum
      )
    }
    out <- lapply(seq_len(N), \(i) {
      cumavg[[i]] |> summarise_with_probs(probs, for_plots) |>
        mutate(
          "{unit}" := treated[i],
          "{time}" := times[seq.int(T_ + 1L, T_total)],
          .before = 1L
        )
    }) |> 
      bind_rows() |> 
      mutate(variable = "Cumulative post-treatment effect")
  }
  if (type == "time") {
    if (average && N > 1) {
      out <- effect |> 
        rvar_apply(1, rvar_mean) |> 
        summarise_with_probs(probs, for_plots) |>
        mutate(
          "{time}" := times,
          .before = 1L
        ) |>
        mutate(variable = "Average treatment effect")
    } else {
      out <- effect |>
        summarise_with_probs(probs, for_plots) |>
        mutate(
          "{unit}" := rep(treated, each = length(times)),
          "{time}" := rep(times, times = length(treated)),
          .before = 1L
        ) |> 
        mutate(variable = "Treatment effect")
    }
  }
  out
}
