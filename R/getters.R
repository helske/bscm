#' Internal getters
#' @noRd
get_treated <- \(x) x$setup$treated
get_donors <- \(x) x$setup$donors
get_treatment <- \(x) x$setup$treatment
get_time <- \(x) x$setup$time
get_unit <- \(x) x$setup$unit
get_times <- \(x) x$setup$times
get_T_pre <- \(x) stats::setNames(x$setup$T_pre, x$setup$treated)
get_T_total <- \(x) x$setup$T_total
get_N <- \(x) length(x$setup$treated)
get_J <- \(x) length(x$setup$donors)
has_intercept <- \(x) x$setup$has_intercept
has_predictors <- \(x) length(x$setup$beta_names) > 0L
has_tv_coefs <- \(x) length(x$setup$gamma_names) > 0L
get_predictors <- \(x) x$setup$predictors

#' Get the value of kappa
#' @param x \[`bscmfit`]\cr The output returned by the [bscm()].
#' @return A value of \eqn{\kappa} used in model estimation.
get_kappa <- function(x) {
  stopifnot_(
    inherits(x, "bscmfit"),
    "Argument {.arg x} should be an object of class {.cls bscmfit}."
  )
  x$setup$kappa
}

#' Extract the name of the outcome variable from formula object
#' @noRd
get_outcome <- function(x) {
  
  if (inherits(x, "bscmfit")) {
    return(x$setup$outcome)
  }
  stopifnot_(
    inherits(x, "formula"), 
    "Argument {.arg formula} must be a {.cls formula} object."
  )
  stopifnot_(
    identical(length(x), 3L), 
    "Argument {.arg formula} must be a {.cls formula} object with an outcome 
    variable on the left-hand side."
  )
  y <- all.vars(x[[2]])
  stopifnot_(
    identical(length(y), 1L), 
    "Argument {.arg formula} must be a {.cls formula} object with one outcome 
    variable on the left-hand side.")
  y
}
