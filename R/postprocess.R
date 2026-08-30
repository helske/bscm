#' Common arguments of the posterior summary methods
#'
#' Documents the arguments shared by the methods which
#' return posterior draws, or posterior summaries, of a quantity computed
#' from a `bscmfit` object.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param summary \[`logical(1)`]\cr If `TRUE` (the default), returns posterior
#'   mean, standard deviation, posterior quantiles (as defined by the
#'   `probs` argument), and MCMC convergence measures.
#'   If `FALSE`, returns the posterior draws instead.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#'   Default is `c(0.025, 0.975)`.
#' @param average \[`logical(1)`]\cr If `TRUE`, the quantity is averaged over
#'   the treated units within each posterior draw in case of multiple treated
#'   units. If `FALSE` (the default), unit-specific values are returned.
#' @param ... Ignored.
#' @name bscm_postprocessing
#' @keywords internal
NULL

#' Replace an rvar column by posterior summaries
#'
#' Summarises the `rvar` column `column` of `d` with
#' `summarise_with_probs()`, and removes the original column together with
#' the `variable` column produced by `posterior::summarise_draws()`.
#'
#' @param d A `tibble` with an `rvar` column named `column`.
#' @param column \[`character(1)`]\cr Name of the `rvar` column.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#' @param for_plots \[`logical(1)`]\cr See `summarise_with_probs()`.
#' @noRd
summarise_column <- function(d, column, probs, for_plots = FALSE) {
  # the summaries are bound rather than mutated in, so that a `variable`
  # column of `d` is not overwritten by the one `summarise_draws()` adds
  summaries <- summarise_with_probs(d[[column]], probs, for_plots) |>
    dplyr::select(-"variable")
  d |>
    dplyr::select(-dplyr::all_of(column)) |>
    dplyr::bind_cols(summaries)
}

#' Extract posterior draws of a single model parameter as an rvar
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param parameter \[`character(1)`]\cr Name of the model parameter.
#' @noRd
rvars_of <- function(x, parameter) {
  as_draws_rvars(x, parameter)[[parameter]]
}

#' Check that an argument is a single logical value
#'
#' @param x Value to check.
#' @param param \[`character(1)`]\cr Name of the argument being checked.
#' @param call The calling environment reported in the error message.
#' @noRd
check_flag <- function(x, param, call = rlang::caller_env()) {
  stopifnot_(
    checkmate::test_flag(x),
    "Argument {.arg {param}} must be a single {.cls logical} value.",
    call = call
  )
}

#' Check that the model fit contains the original data
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param param \[`character(1)`]\cr Name of the argument being checked.
#' @param call The calling environment reported in the error message.
#' @noRd
check_has_data <- function(x, param, call = rlang::caller_env()) {
  stopifnot_(
    !is.null(x$data),
    "The model fit {.arg {param}} does not contain the original data. You
    probably used {.fun bscm} with {.arg save_data = FALSE}?",
    call = call
  )
}

#' Long-format grid of treated units and time points
#'
#' Builds the tibble skeleton shared by the post-processing methods that
#' return one value per treated unit and time point: the unit, time, and
#' treatment indicator columns, in the row order used by the Stan output
#' (all time points of the first treated unit, then the second, and so on).
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param event_time \[`logical(1)`]\cr If `TRUE`, the output also contains
#'   a `time_since_treatment` column, counting time points relative to the
#'   last pre-treatment time point of each treated unit.
#' @noRd
treated_grid <- function(x, event_time = FALSE) {
  T_pre <- get_T_pre(x)
  T_total <- get_T_total(x)
  treated <- get_treated(x)
  unit <- get_unit(x)
  time <- get_time(x)
  treatment <- get_treatment(x)
  d <- dplyr::tibble(
    "{unit}" := rep(treated, each = T_total),
    "{time}" := rep(get_times(x), times = get_N(x)),
    time_since_treatment = unlist(
      lapply(treated, \(i) seq_len(T_total) - T_pre[i] - 1L)
    ),
    "{treatment}" := unlist(
      lapply(treated, \(i) rep(0:1, times = c(T_pre[i], T_total - T_pre[i])))
    )
  )
  if (!event_time) {
    d <- d |> dplyr::select(-"time_since_treatment")
  }
  d
}

#' Check that `probs` defines a single posterior interval
#'
#' Used by the plotting functions, which draw exactly one ribbon and
#' therefore need exactly two probabilities.
#'
#' @param probs \[`numeric(2)`]\cr Probabilities to check.
#' @param call The calling environment reported in the error message.
#' @noRd
interval_probs <- function(probs, call = rlang::caller_env()) {
  stopifnot_(
    checkmate::test_numeric(
      probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      len = 2L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector of length 2 with
    values between 0 and 1.",
    call = call
  )
  sort(probs)
}

#' Check that the model contains predictors
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param call The calling environment reported in the error message.
#' @noRd
check_has_predictors <- function(x, call = rlang::caller_env()) {
  stopifnot_(
    has_predictors(x),
    "The model does not contain any predictors.",
    call = call
  )
}

#' Check that an argument is a single value between 0 and 1
#'
#' @param x Value to check.
#' @param call The calling environment reported in the error message.
#' @noRd
check_alpha <- function(x, call = rlang::caller_env()) {
  stopifnot_(
    checkmate::test_number(x, lower = 0.0, upper = 1.0),
    "Argument {.arg alpha} must be a single {.cls numeric} value between 0
    and 1.",
    call = call
  )
}

#' Resolve and check the `unit` argument
#'
#' Returns the treated units to use: all of them when `unit` is `NULL`.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param unit \[`character()`]\cr Treated units, or `NULL` for all of them.
#' @param call The calling environment reported in the error message.
#' @noRd
check_units <- function(x, unit, call = rlang::caller_env()) {
  treated <- get_treated(x)
  if (is.null(unit)) {
    return(treated)
  }
  stopifnot_(
    all(unit %in% treated),
    "Argument {.arg unit} must be one of the treated units: {.val {treated}}.",
    call = call
  )
  treated[match(unit, treated)]
}

#' Check that `max_lag` is a single positive integer
#'
#' @param max_lag \[`integer(1)`]\cr Value to check.
#' @param call The calling environment reported in the error message.
#' @noRd
check_max_lag <- function(max_lag, call = rlang::caller_env()) {
  stopifnot_(
    checkmate::test_count(max_lag),
    "Argument {.arg max_lag} must be a single positive integer.",
    call = call
  )
}

#' Warn about the deprecated `plot` argument
#'
#' The covariate methods used to draw a plot as a side effect. The argument is
#' absorbed by `...`, so it is detected there.
#'
#' @param ... Arguments passed to the calling method.
#' @param replacement \[`character(1)`]\cr Name of the plotting function that
#'   replaces the removed argument.
#' @noRd
warn_deprecated_plot <- function(..., replacement) {
  if ("plot" %in% ...names()) {
    warning_(c(
      "Argument {.arg plot} is deprecated and no longer has any effect.",
      i = "Use {.fun {replacement}} to plot the results."
    ))
  }
}

#' Convert stored generated quantities to an rvar
#'
#' The generated quantities are stored as a draws by observation matrix, with
#' the draws of each chain stored consecutively. Restoring that chain
#' structure makes the convergence diagnostics of the derived quantities
#' chain-aware. The column names carry the dimensions of the quantity, and are
#' what `posterior::as_draws_rvars()` uses to rebuild them.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param draws A `matrix` of posterior draws stored in the model fit.
#' @noRd
gq_to_rvar <- function(x, draws) {
  stopifnot_(
    identical(nrow(draws), ndraws(x)),
    "The stored posterior draws do not match the number of posterior draws of
    the model fit."
  )
  arr <- array(
    draws,
    dim = c(ndraws(x) / nchains(x), nchains(x), ncol(draws)),
    dimnames = list(NULL, NULL, colnames(draws))
  )
  posterior::as_draws_rvars(arr)[[1L]]
}
