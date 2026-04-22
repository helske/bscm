#' Check arguments of `bscm` call
#'
#' @inheritParams bscm
#' @noRd
check_bscm_arguments <- function(formula, data, treatment, time, unit, 
                                 kappa, r_ess, mcmc_diagnostics, save_data,
                                 priors) {
  
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
    !missing(treatment),
    "Argument {.arg treatment} is missing."
  )
  stopifnot_(
    checkmate::test_string(treatment),
    "Argument {.arg treatment} must be a single character string defining the 
    treatment variable."
  )
  stopifnot_(
    !(treatment %in% all.vars(formula)),
    "Argument {.arg formula} must not contain the treatment variable 
    {.val {treatment}}."
  )
  stopifnot_(
    !is.null(data[[treatment]]),
    "Can't find treatment variable {.var {treatment}} in {.arg data}."
  )
  stopifnot_(
    (checkmate::test_integerish(
      data[[treatment]], lower = 0, upper = 1, any.missing = FALSE
    ) || checkmate::test_logical(
      data[[treatment]], any.missing = FALSE
    )) && length(unique(data[[treatment]])) == 2L,
    "Variable {.arg {treatment}} in {.arg data} should contain either logical 
    or binary values indicating the pre- and post-treatment time points."
  )
  stopifnot_(
    checkmate::test_string(time),
    "Argument {.arg time} must be a single character string defining the 
    time index variable."
  )
  stopifnot_(
    !is.null(data[[time]]),
    "Can't find time index variable {.var {time}} in {.arg data}."
  )
  stopifnot_(
    checkmate::test_numeric(data[[time]], any.missing = FALSE),
    "Time index variable {.arg {time}} must be of type {.cls numeric} and 
    cannot contain missing values."
  )
  stopifnot_(
    checkmate::test_string(unit),
    "Argument {.arg unit} must be a single character string defining the 
    unit index variable."
  )
  stopifnot_(
    !is.null(data[[unit]]),
    "Can't find unit index variable {.var {unit}} in {.arg data}."
  )
  stopifnot_(
    checkmate::test_number(kappa, lower = 0, finite = TRUE, null.ok = TRUE),
    "Argument {.arg kappa} should be single positive number defining the 
      scale parameter of the logistic normal prior of donor weights."
  )
  stopifnot_(
    checkmate::test_number(
      r_ess, lower = 0, upper = 1, finite = TRUE, null.ok = TRUE
    ),
    "Argument {.arg r_ess} should be single number between 0 and 1 defining the 
      prior relative effective number of donors."
  )
  stopifnot_(
    checkmate::test_flag(mcmc_diagnostics),
    "Argument {.arg mcmc_diagnostics} must be a single {.cls logical} value."
  )
  stopifnot_(
    identical(priors, "default"),
    "Argument {.arg priors} is not equal to {.val 'default'}. Only default 
    priors for intercept and sigma are currently supported."
  )
  stopifnot_(
    checkmate::test_flag(save_data),
    "Argument {.arg save_data} must be a single {.cls logical} value."
  )
}

test_probs <- function(probs) {
  stopifnot_(
    checkmate::test_numeric(
      probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      min.len = 0L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values between
     0 and 1."
  )
}