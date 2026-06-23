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
has_intercept <- \(x) x$setup$has_icpt
has_predictors <- \(x) x$setup$has_x
has_tv_coefs <- \(x) x$setup$has_w
get_predictors <- \(x) x$setup$predictors
has_ar1_error <- \(x) x$setup$has_ar1

#' @noRd
get_stan_y <- function(x) {
  unit <- get_unit(x)
  treated <- get_treated(x)
  outcome <- get_outcome(x)
  T_total <- get_T_total(x)
  x$data |>
    dplyr::filter(.data[[unit]] %in% .env$treated) |>
    dplyr::pull(.data[[outcome]]) |>
    matrix(nrow = T_total)
}
#' @noRd
get_stan_Z <- function(x) {
  unit <- get_unit(x)
  donors <- get_donors(x)
  outcome <- get_outcome(x)
  T_total <- get_T_total(x)
  x$data |>
    dplyr::filter(.data[[unit]] %in% .env$donors) |>
    dplyr::pull(.data[[outcome]]) |>
    matrix(nrow = T_total)
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
    variable on the left-hand side."
  )
  y
}
