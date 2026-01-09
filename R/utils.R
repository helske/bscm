#' Helper function for error messages
#' @noRd
stopifnot_ <- function (cond, message, ..., call = rlang::caller_env()) {
  if (!cond) {
    cli::cli_abort(message, ..., .envir = parent.frame(), 
                   call = call)
  }
}
#' Extract the name of the outcome variable from formula object
#' @noRd
get_outcome <- function(x) {
  stopifnot_(
    inherits(x, "formula"), 
    "{.arg formula} must be a {.cls formula} object."
    )
  stopifnot_(
    identical(length(x), 3L), 
    "{.arg formula} must contain the outcome variable on 
    the left-hand side of the {.cls formula} object."
    )
  y <- all.vars(x[[2]])
  stopifnot_(
    length(y) == 1L, 
    "{.arg formula} must be a {.cls formula} object with one outcome variable 
    on the left-hand side.")
  y
}
