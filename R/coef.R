#' Extract regression coefficients of a Bayesian synthetic control model
#'
#' @inheritParams bscm_postprocessing
#' @param type \[`character()`]\cr Type of coefficients to return. Should be
#'   one or more of `"alpha"` (intercepts), `"beta"` (regression coefficients),
#'   `"gamma"` (time-varying regression coefficients), `"kappa"`
#'   (SDs of time-varying coefficients), and
#'   `"rho"` (autoregressive coefficients of the residuals). The default `NULL`
#'   corresponds to all terms above contained in the model.
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
  N <- get_N(object)
  available <- c(
    alpha = has_intercept(object),
    beta = has_predictors(object),
    gamma = has_tv_coefs(object),
    kappa = has_tv_coefs(object),
    rho = has_ar1_error(object)
  )
  available_types <- names(available)[available]
  stopifnot_(
    length(available_types) > 0L,
    "The model does not contain {cli::qty(N)} {?an intercept/intercepts}, any 
    predictors, nor autoregressive coefficients."
  )
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  if (is.null(type)) {
    type <- available_types
  } else {
    type <- try_(match.arg(type, names(available), several.ok = TRUE))
    stopifnot_(
      !inherits(type, "try-error"),
      'Argument {.arg type} must be a subset of
    {{"alpha", "beta", "gamma", "kappa", "rho"}}'
    )

    unavailable <- setdiff(type, available_types)
    stopifnot_(
      length(unavailable) == 0L,
      "The model does not contain the following parameter type{?s}: {unavailable}."
    )
  }
  d <- stats::setNames(
    lapply(type, \(par) coef_draws(object, par)),
    type
  )
  if (summary) {
    d <- stats::setNames(
      lapply(type, \(par) summarise_column(d[[par]], par, probs)),
      type
    )
    # intercepts and fixed effects first, then the time-varying components
    out <- d[intersect(c("alpha", "beta", "rho", "kappa", "gamma"), type)] |>
      dplyr::bind_rows() |>
      dplyr::relocate(
        dplyr::any_of(c("parameter", get_unit(object), get_time(object))),
        .before = 1L
      ) |>
      dplyr::rename(variable = "parameter")
  } else {
    out <- d[intersect(names(available), type)]
  }
  out
}

#' Posterior draws of the regression coefficients of a BSCM
#'
#' Returns the posterior draws of one parameter type as a tibble with a
#' `parameter` column and an rvar column named after the parameter type.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param type \[`character(1)`]\cr One of `"alpha"`, `"beta"`, `"gamma"`,
#'   `"kappa"`, and `"rho"`.
#' @noRd
coef_draws <- function(x, type) {
  unit <- get_unit(x)
  switch(
    type,
    alpha = dplyr::tibble(
      parameter = "Intercept",
      "{unit}" := get_treated(x),
      alpha = rvars_of(x, "alpha")
    ),
    beta = dplyr::tibble(
      parameter = paste0("beta_", x$setup$beta_names),
      beta = rvars_of(x, "beta")
    ),
    gamma = dplyr::tibble(
      parameter = paste0(
        "gamma_",
        rep(x$setup$gamma_names, each = get_T_total(x))
      ),
      "{get_time(x)}" := rep(
        get_times(x),
        times = length(x$setup$gamma_names)
      ),
      gamma = c(rvars_of(x, "gamma"))
    ),
    kappa = dplyr::tibble(
      parameter = paste0("kappa_", x$setup$gamma_names),
      kappa = rvars_of(x, "kappa")
    ),
    rho = dplyr::tibble(
      parameter = "rho",
      "{unit}" := get_treated(x),
      rho = rvars_of(x, "rho")
    )
  )
}
