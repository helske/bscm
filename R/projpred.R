#' Create a projpred reference model from a BSCM fit
#'
#' Creates a `refmodel` object from `bscmfit` to be used with `projpred`
#' package. This enables the usage of [projpred::varsel()] and
#' [projpred::cv_varsel()] for donor selection. Note that for predictions with
#' `newdata`, you need to call [proj_predict_bscm()] which is a wrapper of
#' [projpred::proj_predict()] to which converts the data to the wide format
#' used by `projpred`.
#'
#' @details
#'
#' Projection are based on only on the pretreatment period.
#'
#' This function is experimental and currently only a single treated unit is
#' supported.
#'
#' For projpred integration, donors are treated as separate "predictors" in the
#' formula. Currently only supported model is one without extra predictors.
#' There are also other restrictions, namely lack of support for K-fold.
#'
#' @param object A `bscmfit` object from [bscm()].
#' @param ... Additional arguments passed to [projpred::init_refmodel()].
#' @return An object of class `refmodel`.
#'
#' @examples
#' \donttest{
#' fit <- bscm(
#'   y ~ 1, treatment = "treatment", time = "time", unit = "id",
#'   data = single_treated, refresh = 0, chains = 2, cores = 2
#' )
#' refmodel <- get_refmodel(fit)
#' vs <- projpred::varsel(refmodel)
#' plot(vs, stats = c("elpd", "rmse"), alpha = 0.05)
#' # predictions using the projection:
#' predictions <- projpred::proj_predict(vs)
#' # posterior mean of the predictions
#' colMeans(predictions)
#' }
#' @seealso[projpred::varsel()], [projpred::cv_varsel()]
#' @aliases get_refmodel
#' @export
#' @export get_refmodel
get_refmodel.bscmfit <- function(object, ...) {
  stopifnot_(
    !is.null(object$data),
    "The model fit {.arg object} does not contain the original data. Refit the 
    model with {.arg save_data = FALSE}."
  )

  stopifnot_(
    length(get_predictors(object)) == 0L,
    "Reference model construction currently only supports BSCM without
    additional predictors."
  )
  stopifnot_(
    get_N(object) == 1L,
    "Reference model construction currently only supports BSCM with single 
    treated unit."
  )
  treated <- get_treated(object)
  donors <- original_donor_names <- get_donors(object)
  unit <- get_unit(object)
  time <- get_time(object)
  outcome <- get_outcome(object)

  units <- unique(object$data[[unit]])
  proj_units <- make.names(units, unique = TRUE)

  if (!identical(units, proj_units)) {
    treated_idx <- which(units == treated)
    treated <- proj_units[treated_idx]
    donors <- proj_units[-treated_idx]
  }

  T_pre <- get_T_pre(object)
  proj_data <- as_proj_data(
    object$data,
    unit,
    time,
    outcome,
    units,
    proj_units
  )
  proj_data <- proj_data[seq_len(T_pre), , drop = FALSE]

  proj_formula <- stats::reformulate(
    donors,
    response = treated,
    intercept = has_intercept(object)
  )

  ref_predfun <- function(fit, newdata = NULL) {
    if (!is.null(newdata)) {
      stop(
        paste0(
          "Function `ref_predfun()` for `bscmfit` projection does not ",
          "support `newdata`. This should only be needed for K-fold CV, ",
          "which is not currently supported."
        )
      )
    }
    t(fit$y_mean[, seq_len(T_pre), drop = FALSE])
  }

  proj_predfun <- function(fits, newdata) {
    newdata <- as_proj_data(
      newdata,
      unit,
      time,
      outcome,
      units,
      proj_units
    )
    preds <- lapply(fits, \(fit) {
      if (length(fit$donors) == 0L) {
        return(rep(fit$alpha, nrow(newdata)))
      }
      Z <- as.matrix(newdata[, fit$donors, drop = FALSE])
      fit$alpha + drop(Z %*% fit$omega)
    })
    do.call(cbind, preds)
  }

  extract_model_data <- function(
    object,
    newdata,
    wrhs = NULL,
    orhs = NULL,
    extract_y = TRUE
  ) {
    if (is.null(newdata)) {
      newdata <- object$proj$data
    }
    N <- nrow(newdata)

    list(
      y = if (extract_y) newdata[[treated]],
      weights = rep(1, N),
      offset = rep(0, N)
    )
  }

  # Custom div_minimizer for simplex constraint
  div_minimizer <- function(
    formula,
    data,
    family,
    weights,
    projpred_var,
    projpred_ws_aug,
    verbose_divmin = FALSE,
    ...
  ) {
    donors_prj <- all.vars(formula[-2L])
    responses <- all.vars(formula[[2]])
    # Intercept-only model
    if (length(donors_prj) == 0) {
      return(
        lapply(
          responses,
          \(y) {
            list(
              alpha = mean(data[[y]]),
              omega = numeric(0),
              donors = integer(0)
            )
          }
        )
      )
    }

    Z_prj <- as.matrix(data[, donors_prj, drop = FALSE])
    J <- ncol(Z_prj)
    X <- cbind(1, Z_prj)
    Dmat <- crossprod(X) + diag(1e-8, J + 1)

    Amat <- cbind(
      c(0, rep(1, J)),
      rbind(0, diag(J))
    )
    bvec <- c(1, rep(0, J))

    #Solve QP for one response column
    solve_qp_single <- function(y_target) {
      dvec <- crossprod(X, y_target)
      qp_result <- tryCatch(
        quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1),
        error = function(e) NULL
      )

      if (is.null(qp_result)) {
        alpha_proj <- mean(y_target)
        omega_proj <- rep(1 / J, J)
      } else {
        alpha_proj <- qp_result$solution[1]
        omega_proj <- pmax(qp_result$solution[-1], 0)
        omega_proj <- omega_proj / sum(omega_proj)
      }

      list(
        alpha = alpha_proj,
        omega = omega_proj,
        donors = donors_prj
      )
    }
    lapply(responses, \(y) solve_qp_single(data[[y]]))
  }

  cvrefbuilder <- function(cvfit) {
    stop(
      "K-fold cross-validation is not supported for bscm models. ",
      "Use cv_varsel() with cv_method = 'LOO' instead."
    )
  }
  object$proj <- list(
    data = proj_data,
    formula = proj_formula,
    original_donor_names = stats::setNames(original_donor_names, donors),
    unit = unit,
    time = time,
    outcome = outcome,
    units = units,
    proj_units = proj_units
  )
  projpred::init_refmodel(
    object = object,
    data = proj_data,
    formula = proj_formula,
    family = stats::gaussian(),
    ref_predfun = ref_predfun,
    div_minimizer = div_minimizer,
    proj_predfun = proj_predfun,
    extract_model_data = extract_model_data,
    cvrefbuilder = cvrefbuilder,
    dis = c(as.matrix(object$stanfit, pars = "sigma")),
    ...
  )
}

#' Predictions from a projected BSCM model
#'
#' Convenience wrapper around [projpred::proj_predict()] that accepts data in
#' the original long format used by [bscm()]. If `newdata` is supplied in long
#' format, it is converted internally to the wide format for `projpred`.
#'
#' @param object A projection object or variable selection object returned by
#'   [projpred::project()], [projpred::varsel()], or [projpred::cv_varsel()].
#' @param newdata Optional new data used for predictions. The outcome, unit,
#'   and time variables should have should have names matching the original
#'   data. If `NULL` (the default), predictions are made for the pre-treatment
#'   period of the original data used to fit the model.
#' @param ... Additional arguments passed to [projpred::proj_predict()].
#' @return The output of [projpred::proj_predict()].
#' @export
proj_predict_bscm <- function(object, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    setup <- object$refmodel$fit$proj
    newdata <- as_proj_data(
      newdata,
      unit = setup$unit,
      time = setup$time,
      outcome = setup$outcome,
      units = setup$units,
      proj_units = setup$proj_units
    )
  }

  projpred::proj_predict(
    object = object,
    newdata = newdata,
    ...
  )
}

#' Convert long-format panel data to wide format for projpred
#' @noRd
as_proj_data <- function(
  data,
  unit,
  time,
  outcome,
  units,
  proj_units
) {
  if (!unit %in% names(data)) {
    return(data)
  }
  data[[unit]] <- proj_units[match(data[[unit]], units)]
  data |>
    dplyr::arrange(.data[[time]]) |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(time),
      names_from = dplyr::all_of(unit),
      values_from = dplyr::all_of(outcome)
    ) |>
    dplyr::select(-dplyr::all_of(time))
}
