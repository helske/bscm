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

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

# sum(Q * x) = 0
build_Q <- function(J) {
  k <- seq_len(J - 1)
  Q <- matrix(0, J, J - 1)
  
  for (j in k) {
    Q[seq_len(j), j] <- 1 / sqrt(j * (j + 1))
    Q[j + 1, j] <- -sqrt(j / (j + 1))
  }
  Q
}
