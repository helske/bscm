#' Extract regression coefficients of a Bayesian synthetic control model
#'

#' @param object \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the model coefficients.
#' @aliases coef
#' @export
coef.bscmfit <- function(object, probs = c(0.025, 0.975), ...) {
  test_probs(probs)
  all_pars <- get_stanfit(object)@model_pars
  pars <- intersect(c("alpha", "beta"), all_pars)
  N <- get_N(object)
  stopifnot_(
    length(pars) > 0L,
    "The model does not contain {cli::qty(N)} {?an intercept/intercepts} or any 
    predictors."
  )
  treated <- get_treated(object)
  unit <- get_unit(object)
  d_alpha <- d_beta <- NULL
  if ("alpha" %in% pars) {
    d_alpha <- as_draws(object, "alpha") |> 
      summarise_with_probs(probs = probs) |> 
      mutate("{unit}" := treated, .before = 1L) |> 
      mutate(variable = "Intercept")
  }
  if ("beta" %in% pars) {
    d_beta <- as_draws(object, "beta") |> 
      summarise_with_probs(probs = probs) |> 
      mutate("{unit}" := NA, .before = 1L) |> 
      mutate(variable = paste0("Coef_", object$setup$coef_names))
  } 
  bind_rows(d_alpha, d_beta)
}

