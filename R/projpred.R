#' Create a projpred reference model from a BSCM fit
#'
#' Creates a `refmodel` object from `bscmfit` to be used with `projpred`
#' package. This enables the usage of [projpred::varsel()] and
#' [projpred::cv_varsel()] for donor selection.
#'
#' @details
#'
#' Projection are based on only on the pretreatment period.
#'
#' For projpred integration, donors are treated as separate "predictors" in the
#' formula. Currently only supported model is one without extra predictors.
#' There are also other restrictions, namely lack of support for K-fold. For
#' obtaining predictions based on the projected model, the `newdata` argument
#' needs the data in wide format, see examples.
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
  units <- unique(object$data[[unit]])
  safe_names <- make.names(units, unique = TRUE)
  if (!identical(units, safe_names)) {
    treated_idx <- which(units == treated)
    treated <- safe_names[treated_idx]
    donors <- safe_names[-treated_idx]
    units <- safe_names
    object$setup$donors <- donors
    object$setup$treated <- treated
  }
  T_total <- get_T_total(object)
  T_pre <- get_T_pre(object)
  has_intercept <- has_intercept(object)
  outcome <- get_outcome(object)
  proj_data <- object$data |>
    dplyr::pull(.data[[outcome]]) |>
    matrix(nrow = T_total) |>
    as.data.frame() |>
    stats::setNames(units)
  proj_data <- proj_data[1:T_pre, , drop = FALSE]
  rhs <- paste(donors, collapse = " + ")
  if (has_intercept) {
    proj_formula <- stats::as.formula(paste(treated, "~", rhs))
  } else {
    proj_formula <- stats::as.formula(paste(treated, "~ 0 +", rhs))
  }

  ref_predfun <- function(fit, newdata = NULL) {
    if (is.null(newdata)) {
      newdata <- fit$proj$data[, fit$setup$donors, drop = FALSE]
      eta <- t(as.matrix(fit$stanfit, "y_mean"))
      eta <- eta[1:fit$setup$T_pre, , drop = FALSE]
    } else {
      stop(
        paste0(
          "Function `ref_predfun()` for `bscmfit` projection does not support ",
          "newdata argument. This error should only happen in case of ",
          "K-fold CV or when calling predict(), ",
          "which are not (yet) supported. If you see this error and think ",
          "you need such features, file a issue on Github."
        )
      )
    }
    eta
  }

  proj_predfun <- function(fits, newdata) {
    # Get predictions for each fit
    preds <- lapply(fits, \(fit) {
      if (length(fit$donors) == 0) {
        return(rep(fit$alpha, nrow(newdata)))
      }
      Z <- as.matrix(newdata[, fit$donors, drop = FALSE])
      fit$alpha + c(Z %*% fit$omega)
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

    out <- list(
      weights = rep(1, N),
      offset = rep(0, N)
    )
    if (extract_y) {
      out$y <- newdata[[object$setup$treated]]
    }
    out
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
    S_prj <- length(responses)

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
    original_donor_names = stats::setNames(original_donor_names, donors)
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
