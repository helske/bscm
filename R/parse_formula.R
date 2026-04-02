#' Parse a bscm formula to separate time-constant and time-varying
#' terms
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
        w_terms = character(0)
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
  stopifnot_(
    length(tv_call) == 2L,
    "A {.code tv()} term must have a single formula argument."
  )
  tv_formula <- tryCatch(
    eval(tv_call[[2L]]),
    error = function(e) NULL
  )
  stopifnot_(
    inherits(tv_formula, "formula") && length(tv_formula) == 2L,
    "The argument to {.code tv()} must be a one-sided formula."
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
    w_terms = tv_term_labels
  )
}
