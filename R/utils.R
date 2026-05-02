#' Stop if condition is not true
#'
#' @noRd
stopifnot_ <- function(cond, message, ..., call = rlang::caller_env()) {
  if (!cond) {
    cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
  }
}
#' Warn if condition is not true
#'
#' @noRd
warnifnot_ <- function(cond, message, ...) {
  if (!cond) {
    cli::cli_warn(message, ..., .envir = parent.frame())
  }
}
#' Helper function for warning messages
#'
#' @param message See [cli::cli_warn()].
#' @param ... See [cli::cli_warn()].
#' @noRd
warning_ <- function(message, ...) {
  cli::cli_warn(message, ..., .envir = parent.frame())
}
#' Silent version of `try`
#'
#' Same as [base::try] but with `silent = TRUE`.
#'
#' @param x Expression to evaluate silently.
#' @noRd
try_ <- function(expr) {
  try(expr, silent = TRUE)
}
#' custom summary function for posterior::summarise_draws()
#' @param x draws object
#' @noRd
summarise_with_probs <- function(x, probs, for_plots = FALSE) {
  if (for_plots) {
    summarise_draws(
      x,
      mean,
      if (length(probs) > 1) ~ quantile2(.x, probs = probs)
    )
  } else {
    summarise_draws(
      x,
      mean,
      sd,
      if (length(probs) > 0) ~ quantile2(.x, probs = probs),
      default_convergence_measures()
    )
  }
}

#' Convert posterior draws to a long data frame
#'
#' @param x A `draws_array` object.
#' @param variable Optional variable names. If supplied, must have length 1 or
#'   equal the number of variables in `x`.
#' @noRd
draws_to_long <- function(x, variable = NULL) {
  d <- as_draws_df(x)
  vars <- posterior::variables(d)
  if (!is.null(variable)) {
    stopifnot_(
      length(variable) %in% c(1L, length(vars)),
      "Argument {.arg variable} must have length 1 or match the number of variables."
    )
    variable <- rep(variable, length.out = length(vars))
  } else {
    variable <- vars
  }
  lapply(
    seq_along(vars),
    \(i) {
      var <- vars[i]
      data.frame(
        variable = variable[i],
        value = d[[var]],
        .chain = d$.chain,
        .draw = d$.draw,
        .iteration = d$.iteration
      )
    }
  ) |>
    bind_rows()
}

#' Return posterior output either as summary or long format data frame of draws
#'
#' @param values Posterior values.
#' @param summary Whether to return summaries.
#' @param probs Probabilities for summaries.
#' @param variable Optional variable names.
#' @param for_plots Whether to omit sd/rhat/ess summaries.
#' @noRd
format_posterior_output <- function(
    values,
    summary,
    probs,
    variable = NULL,
    for_plots = FALSE
) {
  if (summary) {
    out <- summarise_with_probs(values, probs, for_plots)
    if (!is.null(variable)) {
      stopifnot_(
        length(variable) %in% c(1L, nrow(out)),
        "Argument {.arg variable} must have length 1 or match the number of rows."
      )
      out <- out |>
        mutate(variable = rep(.env$variable, length.out = nrow(out)))
    }
    return(out)
  }
  draws_to_long(as_draws_array(values), variable = variable)
}

#' Add a unit/time column to summary or draws output
#'
#' @param d A data.frame from ´format_posterior_output()`.
#' @param name Output column name.
#' @param values Column values.
#' @param summary Whether `d` is summary output.
#' @param before Column before which to insert.
#' @param after Column after which to insert.
#' @noRd
add_output_column <- function(
    d,
    name,
    values,
    summary,
    before = 1L
) {
  if (summary) {
    stopifnot_(
      length(values) %in% c(1L, nrow(d)),
      "Length of {.arg values} must be 1 or match output rows."
    )
    vals <- rep(values, length.out = nrow(d))
  } else {
    stopifnot_(
      length(values) > 0L && nrow(d) %% length(values) == 0L,
      "Draw rows are not compatible with {.arg values} length."
    )
    vals <- rep(values, each = nrow(d) / length(values))
  }
  mutate(d, "{name}" := vals, .before = .env$before)
}

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}
