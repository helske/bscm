#' Extract regression coefficients of a Bayesian synthetic control model
#'
#' @inheritParams rmse.bscmfit
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param type \[`character()`]\cr Type of coefficients to return. Should be
#'   one or more of `"alpha"` (intercepts), `"beta"` (regression coefficients),
#'   `"gamma"` (time-varying regression coefficients).
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries (`summary = TRUE`) or
#'   posterior samples (`summary = FALSE`) in long format.
#' @aliases coef
#' @export
#' @examples
#' coef(fit_single_treated)
#' coef(fit_single_treated, type = "beta")
coef.bscmfit <- function(
  object,
  type = c("alpha", "beta", "gamma"),
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  test_summary(summary)
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

  d_alpha <- d_beta <- d_gamma <- d_sigma_gamma <- NULL
  if ("alpha" %in% pars) {
    treated <- get_treated(object)
    unit <- get_unit(object)
    d_alpha <- format_posterior_output(
      as_draws(object, "alpha"),
      summary = summary,
      probs = probs,
      variable = "Intercept"
    ) |>
      add_output_column(
        name = unit,
        values = treated,
        summary = summary
      )
  }
  if ("beta" %in% pars) {
    d_beta <- format_posterior_output(
      as_draws(object, "beta"),
      summary = summary,
      probs = probs,
      variable = paste0("Coef_", object$setup$beta_names)
    )
  }
  if ("gamma" %in% pars) {
    time <- get_time(object)
    times <- get_times(object)
    L <- length(object$setup$gamma_names)
    gamma_names <- object$setup$gamma_names
    gamma_vars <- rep(paste0("gamma_", gamma_names), times = length(times))
    gamma_times <- rep(times, each = L)
    d_gamma <- format_posterior_output(
      as_draws(object, "gamma"),
      summary = summary,
      probs = probs,
      variable = gamma_vars
    )
    d_gamma <- if (summary) {
      add_output_column(
        d_gamma,
        name = time,
        values = gamma_times,
        summary = summary
      )
    } else {
      add_output_column(
        d_gamma,
        name = time,
        values = gamma_times,
        summary = summary
      )
    }
    d_sigma_gamma <- format_posterior_output(
      as_draws(object, "sigma_gamma"),
      summary = summary,
      probs = probs,
      variable = paste0("sigma_gamma_", gamma_names)
    )
  }
  out <- list(
    alpha = d_alpha,
    beta = d_beta,
    gamma = d_gamma,
    sigma_gamma = d_sigma_gamma
  )
  out[lengths(out) > 0]
}
