#' Stop if condition is not true
#' 
#' @noRd
stopifnot_ <- function (cond, message, ..., call = rlang::caller_env()) {
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
summarize_with_probs <- function(x, probs) {
  summarise_draws(
    x,
    mean, sd, 
    ~ quantile2(.x, probs = probs), 
    default_convergence_measures()
  )
}
