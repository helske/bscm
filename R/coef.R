#' Extract regression coefficients of a Bayesian synthetic control model
#'
#' @inheritParams rmse.bscmfit
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param type \[`character()`]\cr Type of coefficients to return. Should be
#'   one or more of `"alpha"` (intercepts), `"beta"` (regression coefficients),
#'   `"gamma"` (time-varying regression coefficients), `"sigma_gamma"` 
#'   (SDs of time-varying coefficients), 
#'   `"rho"` (autoregressive coefficients of the residuals).
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
    type = c("alpha", "beta", "gamma", "sigma_gamma", "rho"),
    summary = TRUE,
    probs = c(0.025, 0.975),
    ...
) {
  test_summary(summary)
  probs <- sort_probs(probs)
  type <- try_(
    match.arg(
      type, c("alpha", "beta", "gamma", "sigma_gamma", "rho"), 
      several.ok = TRUE
    )
  )
  stopifnot_(
    !inherits(type, "try-error"),
    'Argument {.arg type} must a subset of 
    {{"alpha", "beta", "gamma", "sigma_gamma, "rho"}}'
  )
  all_pars <- get_stanfit(object)@model_pars
  pars <- intersect(type, all_pars)
  N <- get_N(object)
  stopifnot_(
    length(pars) > 0L,
    "The model does not contain {cli::qty(N)} {?an intercept/intercepts}, any 
    predictors, nor autoregressive coefficients."
  )
  unit <- get_unit(object)
  time <- get_time(object)
  
  d_alpha <- d_beta <- d_gamma <- d_sigma_gamma <- d_rho <- NULL
  if ("alpha" %in% pars) {
    treated <- get_treated(object)
    d_alpha <- dplyr::tibble(
      parameter = "Intercept",
      "{unit}" := treated,
      alpha = posterior::as_draws_rvars(as_draws(object, "alpha"))$alpha
    )
    if (summary) {
      d_alpha <- d_alpha |>
        dplyr::mutate(summarise_with_probs(.data$alpha, probs)) |>
        dplyr::select(-"alpha", -"variable")
    }
    if (N == 1) d_alpha <- d_alpha |> dplyr::select(-dplyr::all_of(unit))
  }
  if ("beta" %in% pars) {
    d_beta <- dplyr::tibble(
      parameter = paste0("beta_", object$setup$beta_names),
      beta = posterior::as_draws_rvars(as_draws(object, "beta"))$beta
    )
    if (summary) {
      d_beta <- d_beta |>
        dplyr::mutate(summarise_with_probs(.data$beta, probs)) |>
        dplyr::select(-"beta", -"variable")
    }
  }
  if ("gamma" %in% pars) {
    times <- get_times(object)
    gammas <- object$setup$gamma_names
    L <- length(gammas)
    T_total <- get_T_total(object)
    d_gamma <- dplyr::tibble(
      parameter = paste0("gamma_", rep(gammas, each = T_total)),
      "{time}" := rep(times, times = L),
      gamma = c(posterior::as_draws_rvars(as_draws(object, "gamma"))$gamma)
    )
    if (summary) {
      d_gamma <- d_gamma |>
        dplyr::mutate(summarise_with_probs(.data$gamma, probs)) |>
        dplyr::select(-"gamma", -"variable")
    }
  }
  if ("sigma_gamma" %in% pars) {
    d_sigma_gamma <- dplyr::tibble(
      parameter = paste0("sigma_gamma_", object$setup$gamma_names),
      sigma_gamma = posterior::as_draws_rvars(
        as_draws(object, "sigma_gamma")
      )$sigma_gamma
    )
    if (summary) {
      d_sigma_gamma <- d_sigma_gamma |>
        dplyr::mutate(summarise_with_probs(.data$sigma_gamma, probs)) |>
        dplyr::select(-"sigma_gamma", -"variable")
    }
  }
  if ("rho" %in% pars) {
    treated <- get_treated(object)
    d_rho <- dplyr::tibble(
      parameter = "rho",
      "{unit}" := treated,
      rho = posterior::as_draws_rvars(as_draws(object, "rho"))$rho
    )
    if (summary) {
      d_rho <- d_rho |>
        dplyr::mutate(summarise_with_probs(.data$rho, probs)) |>
        dplyr::select(-"rho", -"variable")
    }
    if (N == 1) d_rho <- d_rho |> dplyr::select(-dplyr::all_of(unit))
  }
  dplyr::bind_rows(
    d_alpha,
    d_beta,
    d_rho,
    d_sigma_gamma,
    d_gamma
  ) |>
    dplyr::relocate(dplyr::any_of(c("parameter", unit, time)), .before = 1L) |>
    dplyr::rename(variable = "parameter")
}
