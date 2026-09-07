#' Extract regression coefficients of a Bayesian synthetic control model
#'
#' Returns the coefficients of the model: with `type = "fixed"`, the intercepts
#' and the coefficients of the predictors whose effect is constant in time,
#' and with `type = "varying"`, the total effect of each time-varying
#' predictor at each time point. For a time-varying predictor, the total effect
#' is `beta + gamma_t`, where `beta` is the average level of the effect in the
#' pre-treatment period and `gamma_t` the deviation from it. Every predictor
#' appears in exactly one of the two types.
#'
#' The individual components `beta`, `nu`, `gamma`, `kappa` and `rho` are not
#' returned here; use [as_draws.bscmfit()] to extract them.
#'
#' @inheritParams bscm_postprocessing
#' @param type \[`character(1)`]\cr Type of coefficients to return. Either
#'   `"fixed"` for the intercepts and the time-constant coefficients, or
#'   `"varying"` for the total effect of the time-varying predictors. The
#'   default is `"fixed"`, unless the model contains neither an intercept nor
#'   a predictor with a time-constant coefficient.
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
  has_fixed <- has_intercept(object) || length(fixed_idx(object)) > 0L
  stopifnot_(
    has_fixed || has_tv_coefs(object),
    "The model does not contain {cli::qty(N)} {?an intercept/intercepts} nor
    any predictors."
  )
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  if (is.null(type)) {
    type <- if (has_fixed) "fixed" else "varying"
  } else {
    type <- try_(match.arg(type, c("fixed", "varying")))
    stopifnot_(
      !inherits(type, "try-error"),
      "Argument {.arg type} must be either {.val fixed} or {.val varying}."
    )
    stopifnot_(
      type == "varying" || has_fixed,
      "The model does not contain {cli::qty(N)} {?an intercept/intercepts} nor
      any time-constant coefficients."
    )
    stopifnot_(
      type == "fixed" || has_tv_coefs(object),
      "The model does not contain any time-varying predictors."
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
#' `parameter` column and an `rvar` column.
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param type \[`character(1)`]\cr One of `"fixed"`, `"varying"`, `"alpha"`,
#'   `"beta"`, `"nu"`, `"kappa"`, and `"rho"`.
#' @noRd
coef_draws <- function(x, type) {
  unit <- get_unit(x)
  switch(
    type,
    fixed = {
      idx <- fixed_idx(x)
      dplyr::bind_rows(
        if (has_intercept(x)) {
          dplyr::tibble(
            parameter = "Intercept",
            "{unit}" := get_treated(x),
            coef = rvars_of(x, "alpha")
          )
        },
        if (length(idx) > 0L) {
          dplyr::tibble(
            parameter = paste0("coef_", x$setup$beta_names[idx]),
            coef = rvars_of(x, "beta")[idx]
          )
        }
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
    alpha = dplyr::tibble(
      parameter = "Intercept",
      "{unit}" := get_treated(x),
      alpha = rvars_of(x, "alpha")
    ),
    beta = dplyr::tibble(
      parameter = paste0("beta_", x$setup$beta_names),
      beta = rvars_of(x, "beta")
    ),
    nu = dplyr::tibble(
      parameter = paste0("nu_", x$setup$gamma_names),
      nu = rvars_of(x, "nu")
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

#' Indices of the predictors with a time-constant coefficient
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @noRd
fixed_idx <- function(x) {
  setdiff(seq_along(x$setup$beta_names), x$setup$tv_idx)
}
