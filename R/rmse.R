#' @export
#' @rdname rmse
rmse <- function(x, ...) {
  UseMethod("rmse", x)
}
#' Extract root mean squared errors of a Bayesian synthetic control model
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param average \[`logical(1)`]\cr If `TRUE` (the default), returns the
#' average RMSEs over treated units in case of multiple treated units.
#' If `FALSE`, unit-specific values are returned.
#' @param summary \[`logical(1)`]\cr If `TRUE` (the default), returns posterior
#'   mean, standard deviation, posterior quantiles (as defined by the
#'   `probs` argument), and MCMC convergence measures. 
#'   If `FALSE`, returns the posterior samples instead.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or 
#'   posterior samples (`summary = FALSE`) in long format.
#' @rdname rmse
#' @aliases rmse
#' @export
rmse.bscmfit <- function(
    x, average = TRUE, summary = TRUE, probs = c(0.025, 0.975), ...) {
  test_summary(summary)
  probs <- sort_probs(probs)
  stopifnot_(
    checkmate::test_flag(average),
    "Argument {.arg average} must be a single {.cls logical} value."
  )
  stopifnot_(
    !is.null(x$data),
    "The model fit {.arg x} does not contain the original data. You probably
    used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  N <- get_N(x)
  T_pre <- get_T_pre(x)
  T_total <- get_T_total(x)
  treated <- get_treated(x)
  
  y_rep <- as_draws_rvars(as_draws(x, "y_rep"))$y_rep
  effect <- get_stan_y(x) - y_rep
  
  pre <- post <- ratio <- vector("list", N)
  for (i in seq_len(N)) {
    T_ <- T_pre[treated[i]]
    pre[[i]] <- sqrt(rvar_mean(effect[seq_len(T_), i]^2))
    post[[i]] <- sqrt(rvar_mean(effect[seq.int(T_ + 1L, T_total), i]^2))
    ratio[[i]] <- post[[i]] / pre[[i]]
  }
  if (average && N > 1) {
    vars <- c("Average pre-RMSE", "Average post-RMSE", "Average RMSE ratio")
    values <- lapply(
      list(do.call(c, pre), do.call(c, post), do.call(c, ratio)),
      rvar_mean
    )
    format_posterior_output(
      values,
      summary = summary,
      probs = probs,
      variable = vars
    )
  } else {
    unit <- get_unit(x)
    values <- c(pre, post, ratio)
    vars <- rep(c("Pre-RMSE", "Post-RMSE", "RMSE ratio"), each = N)
    units <- rep(treated, times = 3)
    format_posterior_output(
      values,
      summary = summary,
      probs = probs,
      variable = vars
    ) |>
      add_output_column(
        name = unit,
        values = units,
        summary = summary
      )
  }
}
