#' Stop if Condition is not True
#' 
#' @noRd
stopifnot_ <- function (cond, message, ..., call = rlang::caller_env()) {
  if (!cond) {
    cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
  }
}
#' Warn if Condition is not True
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
#' Silent Version of `try`
#' 
#' Same as [base::try] but with `silent = TRUE`.
#' 
#' @param x Expression to evaluate silently.
#' @noRd
try_ <- function(expr) {
  try(expr, silent = TRUE)
}

