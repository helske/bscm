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
#' The output has the columns `parameter` (`"intercept"` or `"coef"`),
#' `term` (the predictor, not for intercepts), a column named by the unit
#' variable (intercepts) or the time variable (`type = "varying"`), and the
#' posterior summaries or the draws in column `coef`.
#'
#' @inheritParams bscm_postprocessing
#' @param type \[`character(1)`]\cr Type of coefficients to return, either
#'   `"intercept"`, `"fixed"`, or `"varying"`. The default is `"fixed"`, or,
#'   for a model without predictors with a time-constant coefficient, the
#'   first of `"varying"` and `"intercept"` which the model contains.
#' @return A `tibble` of posterior summaries (`summary = TRUE`) or posterior
#'   draws (`summary = FALSE`).
#' @seealso [summary.bscmfit()] for the other model parameters.
#' @aliases coef
#' @export
#' @examples
#' coef(fit_single_treated)
coef.bscmfit <- function(
  object,
  type = NULL,
  summary = TRUE,
  probs = c(0.025, 0.975),
  diagnostics = TRUE,
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
    c(
      "The model contains neither {cli::qty(N)} {?an intercept/intercepts} nor
      any predictors.",
      i = "There are no coefficients to extract."
    )
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
    d <- summarise_column(d, "coef", probs, diagnostics)
  }
  d
}

#' Posterior draws of the coefficients of a BSCM
#'
#' Returns the posterior draws of one type of coefficients returned by
#' [coef.bscmfit()] as a tibble with columns `parameter`, `term` (except for
#' intercepts), a column named by the unit variable (intercepts) or the time
#' variable (`"varying"`), and an `rvar` column named `coef`.
#' 
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param type \[`character(1)`]\cr One of `"intercept"`, `"fixed"`, and
#'   `"varying"`.
#' @noRd
coef_draws <- function(x, type) {
  switch(
    type,
    intercept = parameter_draws(x, "intercept"),
    fixed = parameter_draws(x, "coef") |>
      dplyr::filter(.data$term %in% x$setup$beta_names[fixed_idx(x)]),
    varying = {
      T_total <- get_T_total(x)
      gamma_names <- x$setup$gamma_names
      dplyr::tibble(
        parameter = "coef",
        term = rep(gamma_names, each = T_total),
        "{get_time(x)}" := rep(get_times(x), times = length(gamma_names)),
        # total effect of a time-varying predictor, beta + gamma_t
        value = c(rvars_of(x, "gamma")) +
          rep(rvars_of(x, "beta")[x$setup$tv_idx], each = T_total)
      )
    }
  ) |>
    dplyr::rename(coef = "value")
}

#' Indices of the predictors with a time-constant coefficient
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @noRd
fixed_idx <- function(x) {
  setdiff(seq_along(x$setup$beta_names), x$setup$tv_idx)
}

