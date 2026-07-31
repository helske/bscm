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
#'   functions. Default is `df = 10`.
#' @param noncentered \[`logical(1)`]\cr If `TRUE` (the default), the spline
#'   coefficients are sampled using a non-centered parameterization.
#'   If `FALSE`, a centered parameterization is used.
#'   Depending on the case, one of these might lead to more efficient and
#'   numerically stable sampling, so if you encounter divergences,
#'   try changing this.
#' @param type \[`character(1)`]\cr Prior type for spline coeffients. Either
#'   `"rw1"` or `"rw2"` for first and second order random walks.
#' @return Object of class `tv_term` (a `list`).
#' @export
tv <- function(tv_formula, df = 10, type = "rw1", noncentered = TRUE) {
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
  type <- try_(match.arg(type, c("rw1", "rw2")))
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val rw1} or {.val rw2}."
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
      df = df,
      type = type,
      noncentered = noncentered
    ),
    class = "tv_term"
  )
}

#' Parse a bscm formula to separate time-constant and time-varying terms
#'
#' @param formula A formula potentially containing `tv()`.
#' @return A list with components:
#'  * full_formula`: formula with all terms (for X matrix)
#'  * `tv_formula`: formula with only tv terms (for column name matching), or
#'  `NULL`
#'  * `has_tv`: logical
#'  * `tv_terms`: character vector of time-varying term labels
#'  * `has_intercept`: logical
#'  * `predictors`: character vector of predictor variable names
#'  * `spline_df`: degrees of freedom
#'  * `spline_type`: character indicating the type of spline basis to use
#'  * `noncentered`: logical indicating whether to use non-centered
#'      parameterization for spline coefficients
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
        spline_df = 0,
        spline_type = "ts",
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

  spline_df <- tv_object$df
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
    spline_type = tv_object$type,
    noncentered_xi = tv_object$noncentered
  )
}
