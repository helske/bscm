#' Check Arguments of `bscm` Call
#'
#' @inheritParams bscm
#' @noRd
check_bscm_arguments <- function(formula, data, treatment, time, unit,
                                 priors, intercept, effective_donors, probs) {
  
  stopifnot_(
    !missing(formula),
    "Argument {.arg formula} is missing."
  )
  stopifnot_(
    inherits(formula, "formula"),
    "Argument {.arg formula} must be a {.cls formula} object."
  )
  stopifnot_(
    !missing(data),
    "Argument {.arg data} is missing."
  )
  stopifnot_(
    is.data.frame(data),
    "Argument {.arg data} must be a {.cls data.frame} object."
  )
  stopifnot_(
    checkmate::test_string(x = treatment),
    "Argument {.arg treatment} must be a single character string defining the 
    treatment variable."
  )
  stopifnot_(
    !is.null(data[[treatment]]),
    "Can't find treatment variable {.var {treatment}} in {.arg data}."
  )
  stopifnot_(
    checkmate::test_string(x = time),
    "Argument {.arg time} must be a single character string defining the 
    time index variable."
  )
  stopifnot_(
    !is.null(data[[time]]),
    "Can't find time index variable {.var {time}} in {.arg data}."
  )
  stopifnot_(
    checkmate::test_numeric(data[[time]], any.missing = FALSE),
    "Time index variable {.arg {time}} must be of of type {.cls numeric} and 
    cannot contain missing values."
  )
  stopifnot_(
    checkmate::test_string(x = unit),
    "Argument {.arg unit} must be a single character string defining the 
    unit index variable."
  )
  stopifnot_(
    !is.null(data[[unit]]),
    "Can't find unit index variable {.var {unit}} in {.arg data}."
  )
  stopifnot_(
    identical(priors, "default"),
    "Argument {.arg priors} is not equal to {.val 'default'}. Only default 
    priors for intercept and sigma are currently supported."
  )
  stopifnot_(
    checkmate::test_flag(x = intercept),
    "Argument {.arg intercept} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_integerish(
      x = effective_donors, len = 1, lower = 2, null.ok = TRUE),
    "Argument {.arg effective_donors} must be a single integer defining the 
    prior expected number of effective donors."
  )
  stopifnot_(
    checkmate::test_numeric(
      x = probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values between
     0 and 1."
  )
}