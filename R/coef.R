#' Extract regression coefficients of a Bayesian synthetic control model
#'
#' Returns one type of coefficients of the model: the intercepts
#' (`type = "intercept"`), the coefficients of the predictors whose effect is
#' constant in time (`type = "fixed"`), or the total effect of each
#' time-varying predictor at each time point (`type = "varying"`). For a
#' time-varying predictor, the total effect is `beta + gamma_t`, where `beta`
#' is the average level of the effect in the pre-treatment period and
#' `gamma_t` the deviation from it. Every predictor appears in exactly one of
#' the types `"fixed"` and `"varying"`.
#'
#' The individual components `beta`, `nu`, `gamma`, `kappa` and `rho` are not
#' returned here; use [as_draws.bscmfit()] to extract them.
#'
#' @inheritParams bscm_postprocessing
#' @param type \[`character(1)`]\cr Type of coefficients to return, either
#'   `"intercept"`, `"fixed"`, or `"varying"`. The default is `"fixed"`, or,
#'   for a model without predictors with a time-constant coefficient, the
#'   first of `"varying"` and `"intercept"` which the model contains.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or posterior
#'   draws (`summary = FALSE`).
#' @seealso [as_draws.bscmfit()] for the individual model parameters.
#' @aliases coef
#' @export
#' @examples
#' coef(fit_single_treated)
coef.bscmfit <- function(
  object,
  type = NULL,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  N <- get_N(object)
  # the names are in the order of preference used by the default type
  available <- c(
    fixed = length(fixed_idx(object)) > 0L,
    varying = has_tv_coefs(object),
    intercept = has_intercept(object)
  )
  stopifnot_(
    any(available),
    "The model does not contain {cli::qty(N)} {?an intercept/intercepts} nor
    any predictors."
  )
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  if (is.null(type)) {
    type <- names(available)[available][1L]
  } else {
    type <- try_(match.arg(type, c("intercept", "fixed", "varying")))
    stopifnot_(
      !inherits(type, "try-error"),
      "Argument {.arg type} must be one of {.val intercept}, {.val fixed}, or
      {.val varying}."
    )
    stopifnot_(
      available[type],
      switch(
        type,
        intercept = "The model does not contain
        {cli::qty(N)} {?an intercept/intercepts}.",
        fixed = "The model does not contain any time-constant coefficients.",
        varying = "The model does not contain any time-varying predictors."
      )
    )
  }
  d <- coef_draws(object, type)
  if (summary) {
    d <- d |>
      summarise_column("coef", probs) |>
      dplyr::rename(variable = "parameter")
  }
  d
}

#' Posterior draws of the coefficients of a BSCM
#'
#' Returns the posterior draws of one type of coefficients as a tibble with a
#' `parameter` column and an `rvar` column named `coef`. The types
#' `"intercept"`, `"fixed"`, and `"varying"` are the views returned by
#' [coef.bscmfit()], the rest are the individual model parameters shown by
#' the print method.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param type \[`character(1)`]\cr One of `"intercept"`, `"fixed"`,
#'   `"varying"`, `"beta"`, `"nu"`, `"kappa"`, and `"rho"`.
#' @noRd
coef_draws <- function(x, type) {
  unit <- get_unit(x)
  switch(
    type,
    intercept = dplyr::tibble(
      parameter = "Intercept",
      "{unit}" := get_treated(x),
      coef = rvars_of(x, "alpha")
    ),
    fixed = {
      idx <- fixed_idx(x)
      dplyr::tibble(
        parameter = paste0("coef_", x$setup$beta_names[idx]),
        coef = rvars_of(x, "beta")[idx]
      )
    },
    varying = {
      T_total <- get_T_total(x)
      gamma_names <- x$setup$gamma_names
      dplyr::tibble(
        parameter = paste0("coef_", rep(gamma_names, each = T_total)),
        "{get_time(x)}" := rep(get_times(x), times = length(gamma_names)),
        # total effect of a time-varying predictor, beta + gamma_t
        coef = c(rvars_of(x, "gamma")) +
          rep(rvars_of(x, "beta")[x$setup$tv_idx], each = T_total)
      )
    },
    beta = dplyr::tibble(
      parameter = paste0("beta_", x$setup$beta_names),
      coef = rvars_of(x, "beta")
    ),
    nu = dplyr::tibble(
      parameter = paste0("nu_", x$setup$gamma_names),
      coef = rvars_of(x, "nu")
    ),
    kappa = dplyr::tibble(
      parameter = paste0("kappa_", x$setup$gamma_names),
      coef = rvars_of(x, "kappa")
    ),
    rho = dplyr::tibble(
      parameter = "rho",
      "{unit}" := get_treated(x),
      coef = rvars_of(x, "rho")
    )
  )
}

#' Indices of the predictors with a time-constant coefficient
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @noRd
fixed_idx <- function(x) {
  setdiff(seq_along(x$setup$beta_names), x$setup$tv_idx)
}
