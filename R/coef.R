#' Extract regression coefficients of a Bayesian synthetic control model
#'
#' @inheritParams rmse.bscmfit
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param type \[`character()`]\cr Type of coefficients to return. Should be
#'   one or more of `"alpha"` (intercepts), `"beta"` (regression coefficients),
#'   `"gamma"` (time-varying regression coefficients), `"sigma_gamma"`
#'   (SDs of time-varying coefficients), and
#'   `"rho"` (autoregressive coefficients of the residuals). The default `NULL`
#'   corresponds to all terms above contained in the model.
#' @param ... Ignored.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or
#'   list of `tibbles` of posterior draws (`summary = FALSE`).
#' @aliases coef
#' @export
#' @examples
#' coef(fit_single_treated)
#' coef(fit_single_treated, type = "beta")
coef.bscmfit <- function(
  object,
  type = NULL,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  stopifnot_(
    has_intercept(object) || has_predictors(object) || has_ar1_error(object),
    "The model does not contain {cli::qty(N)} {?an intercept/intercepts}, any 
    predictors, nor autoregressive coefficients."
  )
  test_summary(summary)
  probs <- sort_probs(probs)
  available_types <- intersect(
    c("alpha", "beta", "gamma", "sigma_gamma", "rho"),
    setdiff(
      get_stanfit(object)@model_pars,
      object$setup$excluded_pars
    )
  )
  if (is.null(type)) {
    type <- available_types
  } else {
    type <- try_(
      match.arg(
        type,
        c("alpha", "beta", "gamma", "sigma_gamma", "rho"),
        several.ok = TRUE
      )
    )
    stopifnot_(
      !inherits(type, "try-error"),
      'Argument {.arg type} must be a subset of
    {{"alpha", "beta", "gamma", "sigma_gamma", "rho"}}'
    )

    unavailable <- setdiff(type, available_types)
    stopifnot_(
      length(unavailable) == 0L,
      "The model does not contain the following parameter type{?s}: {unavailable}."
    )
  }
  N <- get_N(object)
  unit <- get_unit(object)
  time <- get_time(object)

  d_alpha <- d_beta <- d_gamma <- d_sigma_gamma <- d_rho <- NULL
  if ("alpha" %in% type) {
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
  if ("beta" %in% type) {
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
  if ("gamma" %in% type) {
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
  if ("sigma_gamma" %in% type) {
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
  if ("rho" %in% type) {
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
  if (summary) {
    out <- dplyr::bind_rows(
      d_alpha,
      d_beta,
      d_rho,
      d_sigma_gamma,
      d_gamma
    ) |>
      dplyr::relocate(
        dplyr::any_of(c("parameter", unit, time)),
        .before = 1L
      ) |>
      dplyr::rename(variable = "parameter")
  } else {
    out <- list(
      alpha = d_alpha,
      beta = d_beta,
      gamma = d_gamma,
      sigma_gamma = d_sigma_gamma
    )
    out <- out[lengths(out) > 0]
  }
  out
}
