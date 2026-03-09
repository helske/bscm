#' Extract regression coefficients of a Bayesian synthetic control model
#'
#' @export
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the coefficients.
coef.bscmfit <- function(object, probs = c(0.025, 0.975), ...) {
  stopifnot_(
    checkmate::test_numeric(
      probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values between
     0 and 1."
  )
  all_pars <- object$stanfit@model_pars
  pars <- intersect(c("alpha", "beta"), all_pars)
  stopifnot_(
    length(pars) > 0L,
    "The model does not contain an intercept or any predictors."
  )
  d_alpha <- d_beta <- NULL
  if ("alpha" %in% pars) {
    d_alpha <- as_draws(object, "alpha") |> 
      summarise_draws(
        mean, sd, 
        ~ quantile2(.x, probs = probs), 
        default_convergence_measures()
      ) |> 
      mutate(variable = "Intercept")
  }
  if ("beta" %in% pars) {
    d_beta <- as_draws(object, "beta") |> 
      summarise_draws(
        mean, sd, 
        ~ quantile2(.x, probs = probs), 
        default_convergence_measures()
      ) |> 
      mutate(variable = paste0("Coef_", object$setup$coef_names))
  } 
  bind_rows(d_alpha, d_beta)
}
