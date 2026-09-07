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
#' @param probs quantile probabilities
#' @param diagnostics if `TRUE` (the default), the summaries also contain the
#'   MCMC convergence measures and the Monte Carlo standard error of the mean
#' @noRd
summarise_with_probs <- function(x, probs, diagnostics = TRUE) {
  posterior::summarise_draws(
    x,
    mean,
    sd,
    if (length(probs) > 0) ~ posterior::quantile2(.x, probs = probs),
    if (diagnostics) posterior::default_convergence_measures(),
    if (diagnostics) "mcse_mean"
  )
}

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}
