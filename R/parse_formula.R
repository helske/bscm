#' Define the time-varying coefficients for BSCM
#'
#' Function `tv()` validates its arguments and returns the necessary
#' variables needed for the construct the time-varying coefficients in
#' [bscm()]. You should not call `tv()` separately, instead it should be part
#' of a formula, e.g., `y ~ x + tv(~ z * w, df = 20)`.
#'
#' @param tv_formula \[`formula`]\cr One-sided formula, e.g., `~ x + z`
#'   defining the the time-varying part of the BSCM formula.
#' @param df \[`integer(1)`]\cr Integer defining the number of spline basis
#'   functions. Default `df = 10` somewhat is arbitrary, so it is good to
#'   test various choices of `df`. Typically it is better to choose too large
#'   than too small value, as the random walk prior of the spline coefficients
#'   will regularize overfitting.
#' @param noncentered \[`logical(1)`]\cr If `TRUE` (the default), the spline
#'   coefficients are sampled using a non-centered parameterization.
#'   If `FALSE`, a centered parameterization is used.
#'   Depending on the case, one of these might lead to more efficient and
#'   numerically stable sampling, so if you encounter divergences,
#'   try changing this.
#' @return Object of class `tv_term` (a `list`).
#' @export
tv <- function(tv_formula, df = 10, noncentered = TRUE) {
  stopifnot_(
    inherits(tv_formula, "formula") && length(tv_formula) == 2L,
    "Argument {.arg tv_formula} of {.code tv()} must be a one-sided formula."
  )
  stopifnot_(
    checkmate::test_int(df, lower = 4),
    "Argument {.arg df} of {.code tv()} must be an integer larger than 4."
  )
  stopifnot_(
    checkmate::test_flag(noncentered),
    "Argument {.arg noncentered} of {.code tv()} must be a logical value."
  )
  tv_terms_obj <- stats::terms(tv_formula, specials = "tv")
  nested_tv <- attr(tv_terms_obj, "specials")$tv
  stopifnot_(
    is.null(nested_tv),
    "Nested {.code tv()} terms are not supported."
  )
  tv_term_labels <- attr(tv_terms_obj, "term.labels")
  stopifnot_(
    length(tv_term_labels) > 0L,
    "The formula inside {.code tv()} must contain at least one term."
  )
  structure(
    list(
      terms = tv_term_labels,
      spline_df = df,
      noncentered_xi = noncentered
    ),
    class = "tv_term"
  )
}

#' Parse a bscm formula to separate time-constant and time-varying terms
#'
#' @param formula A formula potentially containing `tv()`.
#' @return A list with components:
#'   - `full_formula`: formula with all terms (for X matrix)
#'   - `tv_formula`: formula with only tv terms (for column name
#'     matching), or NULL
#'   - `has_tv`: logical
#'   - `tv_terms`: character vector of time-varying term labels
#'   - `has_intercept`: logical
#'   - `predictors`: character vector of predictor variable names
#'   - `spline_df`: degrees of freedom variable passed on to [splines::bs()].
#'   - `noncentered_xi`: logical indicating whether to use non-centered
#'      parameterization for spline coefficients.
#' @noRd
parse_bscm_formula <- function(formula) {
  xt <- stats::terms(formula, specials = "tv")
  tv_idx <- attr(xt, "specials")$tv
  icpt <- as.logical(attr(xt, "intercept"))

  if (is.null(tv_idx)) {
    predictors <- all.vars(formula[-2L])
    return(
      list(
        icpt = icpt,
        predictors = predictors,
        x_formula = formula,
        w_formula = NULL,
        w_terms = character(0),
        df = 0,
        noncentered_xi = TRUE
      )
    )
  }

  stopifnot_(
    length(tv_idx) == 1L,
    "Multiple {.code tv()} terms are not supported. Combine all
    time-varying terms into a single {.code tv()} call."
  )

  xt_vars <- attr(xt, "variables")
  xt_terms <- attr(xt, "term.labels")

  tv_call <- xt_vars[[tv_idx + 1L]]
  tv_object <- eval(tv_call)

  spline_df <- tv_object$spline_df
  tv_term_labels <- tv_object$terms

  # Drop the tv() term from the original formula
  tv_deparse <- deparse1(xt_vars[[tv_idx + 1L]])
  tv_drop_idx <- which(xt_terms == tv_deparse)

  main_xt <- stats::drop.terms(xt, dropx = tv_drop_idx)
  main_terms <- attr(stats::terms(main_xt), "term.labels")
  response <- xt_vars[[2L]]

  # union of main and tv terms
  all_terms <- union(main_terms, tv_term_labels)
  full_formula <- stats::reformulate(
    termlabels = all_terms,
    response = response,
    intercept = icpt
  )

  # tv terms only
  tv_formula <- stats::reformulate(
    termlabels = tv_term_labels,
    response = response,
    intercept = icpt
  )

  predictors <- all.vars(full_formula[-2L])

  list(
    icpt = icpt,
    predictors = predictors,
    x_formula = full_formula,
    w_formula = tv_formula,
    w_terms = tv_term_labels,
    spline_df = spline_df,
    noncentered_xi = tv_object$noncentered_xi
  )
}
