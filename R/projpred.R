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
#' \dontrun{
#' fit <- bscm(
#'   y ~ 1, treatment = "treatment", time = "time", unit = "id",
#'   data = simulated_data, refresh = 0, chains = 2, cores = 2
#' )
#' refmodel <- get_refmodel(fit)
#' vs <- projpred::varsel(refmodel)
#' plot(vs)
#' # predict with new data using the projection:
#' # (here original data for all time points)
#' # simulated_data is ordered by id and time
#' newdata <- matrix(simulated_data$y, nrow = 50)
#' colnames(newdata) <- unique(simulated_data$id)
#' newdata <- as.data.frame(newdata)
#' # returns posterior samples of predictions
#' predictions <- projpred::proj_predict(
#'   vs, newdata = newdata
#' )
#' # posterior mean of the predictions
#' colMeans(predictions)
#' } 
#' @seealso[projpred::varsel()], [projpred::cv_varsel()]
#' @aliases get_refmodel
#' @export get_refmodel
#' @export
get_refmodel.bscmfit <- function(object, ...) {
  
  stopifnot_(
    !is.null(object$data),
    "The model fit does not contain the original data. Refit with 
    {.arg save_data = TRUE}."
  )
  stopifnot_(
    length(object$setup$predictors) == 0,
    "Reference model construction currently only supports BSCM without
    additional predictors."
  )
  unit <- object$setup$unit
  T_total <- object$setup$T_total
  T_pre <- object$setup$T_pre
  has_intercept <- object$setup$has_intercept
  outcome <- object$setup$outcome
  treated <- object$setup$treated
  donors <- object$setup$donors
  
  units <- unique(object$data[[unit]])
  proj_data <- object$data |> 
    pull(.data[[outcome]]) |> 
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
      eta <- t(as.matrix(fit$stanfit, "synthetic_mean"))
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
    preds <- lapply(fits, function(fit) {
      if (length(fit$donors) == 0) {
        return(rep(fit$alpha, nrow(newdata)))
      }
      Z <- as.matrix(newdata[, fit$donors, drop = FALSE])
      fit$alpha + c(Z %*% fit$omega)
    })
    do.call(cbind, preds)
  }
  
  extract_model_data <- function(object, newdata, wrhs = NULL, orhs = NULL,
                                 extract_y = TRUE) {
    
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
  div_minimizer <- function(formula, data, family, weights, 
                            projpred_var, projpred_ws_aug, 
                            verbose_divmin = FALSE, ...) {
    donors_prj <- all.vars(formula[-2L])
    responses <- all.vars(formula[[2]])
    S_prj <- length(responses)
    
    # Intercept-only model
    if (length(donors_prj) == 0) {
      return(
        lapply(
          responses, \(y) list(
            alpha = mean(data[[y]]), omega = numeric(0), donors = integer(0)
          )
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
    stop("K-fold cross-validation is not supported for bscm models. ",
         "Use cv_varsel() with cv_method = 'LOO' instead.")
  }
  object$proj <- list(
    data = proj_data, formula = proj_formula
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
#' Extract Donor Ranking from projpred varsel Object
#' 
#' Convenience function to extract the donor ranking from a `vsel` object
#' returned by [projpred::varsel()] or [projpred::cv_varsel()] when applied to 
#' a `bscmfit` object.
#' 
#' @param x \[`vsel`\] Output from [projpred::varsel()] or 
#' [projpred::cv_varsel()].
#' @return Character vector of donor IDs in order of selection.
#' @export
donor_ranking <- function(x) {
  stopifnot_(
    inherits(x, "vsel"),
    "{.arg x} must be a 'vsel' object from {.fn projpred::varsel} or 
    {.fn projpred::cv_varsel}."
  )
  projpred::ranking(x)$fulldata
}

