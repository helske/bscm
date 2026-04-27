#' Extract regression coefficients of a Bayesian synthetic control model
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param type \[`character()`]\cr Type of coefficients to return. Should be 
#' one or more of `"alpha"` (intercepts), `"beta"` (regression coefficients),
#' `"gamma"` (time-varying regression coefficients).
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the model coefficients.
#' @aliases coef
#' @export
coef.bscmfit <- function(
    object, type = c("alpha", "beta", "gamma"), probs = c(0.025, 0.975), ...) {
  probs <- sort_probs(probs)
  type <- try_(
    match.arg(type, c("alpha", "beta", "gamma"), several.ok = TRUE)
    )
  stopifnot_(
    !inherits(type, "try-error"),
    'Argument {.arg type} must a subset of {{"alpha", "beta", "gamma"}}'
  )
  all_pars <- get_stanfit(object)@model_pars
  pars <- intersect(type, all_pars)
  N <- get_N(object)
  stopifnot_(
    length(pars) > 0L,
    "The model does not contain {cli::qty(N)} {?an intercept/intercepts} or any 
    predictors."
  )
  
 
  d_alpha <- d_beta <- d_gamma <- d_kappa <- NULL
  if ("alpha" %in% pars) {
    treated <- get_treated(object)
    unit <- get_unit(object)
    d_alpha <- as_draws(object, "alpha") |> 
      summarise_with_probs(probs = probs) |> 
      mutate("{unit}" := treated, .before = 1L) |> 
      mutate(variable = "Intercept")
  }
  if ("beta" %in% pars) {
    d_beta <- as_draws(object, "beta") |> 
      summarise_with_probs(probs = probs) |> 
      mutate(variable = paste0("Coef_", object$setup$beta_names))
  }
  if ("gamma" %in% pars) {
    time <- get_time(object)
    times <- get_times(object)
    L <- length(object$setup$gamma_names)
    d_gamma <- as_draws(object, "gamma") |> 
      summarise_with_probs(probs = probs) |> 
      mutate("{time}" := rep(times, each = L), .before = 1L) |> 
      mutate(variable = paste0("gamma_", object$setup$gamma_names))
    d_kappa <- as_draws(object, "kappa") |> 
      summarise_with_probs(probs = probs) |> 
      mutate(variable = paste0("kappa_", object$setup$gamma_names))
  } 
  out <- list(alpha = d_alpha, beta = d_beta, gamma = d_gamma, kappa = d_kappa)
  out[lengths(out) > 0]
}

