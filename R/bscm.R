#' Estimate Bayesian synthetic control model
#'
#' Function `bscm` estimates a Bayesian synthetic control models of varying
#' types.
#'
#' @details
#'
#' To define `formula` in case of no predictors, use
#' `outcome ~ 1` or `outcome ~ 0`. In the former case, as well by the
#' default when using predictors (e.g., `outcome ~ x + z`), the model
#' includes intercept term for each treated unit. Intercept can be omitted by
#' using `0` in the RHS, e.g., `outcome ~ 0` or `outcome ~ 0 + x + z`
#' (equivalently, you can use `-1` in place of `0`). Formula should not contain
#' the variable defining the treatment, which is defined separately using the
#' argument `treatment`. In case the variable is present also in the
#' formula, it is automatically removed.
#'
#' To specify predictors with time-varying coefficients in `formula`, wrap them
#' in `tv()`, e.g., `outcome ~ x + tv(~ z, df = 10, type = "rw1")` defines a
#' model where `x` has a time-constant coefficient and `z` has a time-varying
#' coefficient following penalized cubic spline with `10` spline basis
#' functions and random walk prior on the spline coefficients.
#' Terms inside `tv()` are automatically also included in the
#' time-constant part of the model, since the time-varying coefficients
#' are defined to have zero mean in the pre-treatment period (in case of
#' multiple treated units, minimum period). Note that a common time-varying
#' intercept is omitted as it would cancel out in the linear predictor.
#'
#' Both the time-constant and time-varying regression part is assumed to apply
#' for all treated and donor units with common coefficients. If you
#' want to apply a covariate only for treated units, just set the covariate
#' values to zero for donors. This approach can be also used to vary behaviour
#' of intercept: Defining formula such as `y ~ 0 + intercept + x`, where
#' `intercept` is a name of constant column in the data, will define common
#' intercept for all units, instead of unit-specific intercepts
#' (e.g., fixed effects).
#'
#' When model contains covariates \eqn{X}, their effect is subtracted from
#' donors, i.e., for treated unit \eqn{i},
#' \eqn{y_i \sim N(\alpha_i + X_i\beta + Z^\ast\omega_i, \sigma_i^2)},
#' where column \eqn{j} of \eqn{Z^\ast} is
#' \eqn{z_j - (\alpha_j + X_j \beta)}, i.e., the original donor outcome
#' vector \eqn{z_j} minus the donor-specific intercept \eqn{\alpha_j} and
#' the effect of covariates \eqn{X_j \beta} for that donor. Note that
#' \eqn{\beta} is common across donors and the treated units.
#'
#' Model can contain missing values in the outcome variable of the treated,
#' but not of the donors, nor in the covariates. Missing outcomes are
#' automatically imputed during MCMC sampling under a missing at random
#' (MAR) assumption.
#'
#' @param formula \[`formula`]\cr The model formula containing the outcome
#' variable on the left-hand side and optional time-varying predictors on
#' the right-hand side (RHS) of `~`. See details.
#' @param data \[`data.frame` or an object coercible to one]\cr
#'   The long format data that contains the model variables.
#' @param treatment  \[`character(1)`]\cr Name of the treatment indicator
#'   variable in `data`. Default is `"treatment"`.
#' @param time \[`character(1)`]\cr Name of the time index variable in
#'   `data`. Default is `"time"`.
#' @param unit \[`character(1)`]\cr Name of the variable in `data`
#'   identifying different units. Default is `"id"`.
#' @param error \[`character(1)`]\cr Assumed error structure of the model.
#'   Either `"iid"` (independent and identically distributed errors) or
#'   `"ar1"` (first-order autoregressive process). Default is `"iid"`.
#' @param priors \[`list()`]\cr List of prior definitions for the model
#'   parameters. See [bscm_prior()] for details on the available prior
#' families and how to define them. If `NULL` (the default), default
#' priors are used.
#' @param prior_only \[`logical(1)`]\cr If `TRUE`, samples from prior
#'   predictive distribution. Default is `FALSE`.
#' @param mcmc_diagnostics \[`logical(1)`]\cr If `TRUE` (the default), the
#'   output of [bscm()] includes the results of MCMC diagnostics checks
#'   performed by [check_mcmc_diagnostics.bscmfit()]. Note that regardless
#'   of the value of `mcmc_diagnostics`, [rstan::sampling()] can still
#'   generate warnings regarding convergence and other sampling issues.
#' @param save_data \[`logical(1)`]\cr If `TRUE` (the default), `bscmfit`
#'   object returned by  [bscm()] includes the input `data.frame`
#'   (argument `data`), after dropping unused factor levels and potentially
#'   rearranging data by `unit` and `time` variables.
#' @param compute_predictions \[`logical(1)`]\cr If `TRUE` (the default),
#'   posterior predictive draws (`y_rep`) are computed during sampling. Set to
#'   `FALSE` to skip this computation and reduce memory usage; this is used
#'   internally by [lfo()] when refitting the model repeatedly. Note that
#'   setting this to `FALSE` will cause [treatment_effect()],
#'   [synthetic_control()], [posterior_predict()], [rmse()], [summary()], and
#'   other methods that rely on posterior predictions to fail, so you rarely
#'   want to set this to `FALSE`.
#' @param ... Additional parameters passed on to [rstan::sampling()] to
#'   adjust the sampling options, for example `iter` and `chains`.  Note that
#'   defaults `iter = 5000` and `warmup = 2500` differ from the defaults of
#'   [rstan::sampling()] (which are 2000 and 1000 respectively). Many 
#'   control arguments for Stan can be passed using a named list `control`,
#'   such as `control = list(adapt_delta = 0.95)` which corresponds to the
#'   default `adapt_delta` value of `bscm` when model contains time-varying 
#'   coefficients (otherwise the default is `0.8` as in `rstan`).
#' @return An object of class `bscmfit`.
#' @export
#' @seealso [summary.bscmfit()], [as_draws.bscmfit()], [rstan::sampling()].
#' @examples
#' # skip diagnostics and use small number of iterations for CRAN checks
#' fit <- bscm(
#'   y ~ 1, single_treated, "treatment", "time", "id",
#'   priors = list(omega = dirichlet_pr(0.5)),
#'   chains = 1, cores = 1, refresh = 0, iter = 1000,
#'   mcmc_diagnostics = FALSE
#' )
#' fit
bscm <- function(
  formula,
  data,
  treatment = "treatment",
  time = "time",
  unit = "id",
  error = "iid",
  priors = NULL,
  prior_only = FALSE,
  mcmc_diagnostics = TRUE,
  save_data = TRUE,
  compute_predictions = TRUE,
  ...
) {
  check_bscm_arguments(
    formula,
    data,
    treatment,
    time,
    unit,
    mcmc_diagnostics,
    save_data,
    compute_predictions,
    prior_only
  )
  outcome <- get_outcome(formula)
  parsed_formula <- parse_bscm_formula(formula)
  formula <- parsed_formula$x_formula
  has_icpt <- parsed_formula$icpt
  predictors <- parsed_formula$predictors
  has_x <- length(predictors) > 0
  has_w <- length(parsed_formula$w_terms) > 0
  spline_df <- parsed_formula$spline_df
  knot_spacing <- parsed_formula$knot_spacing
  spline_type <- parsed_formula$spline_type
  noncentered_xi <- parsed_formula$noncentered_xi
  stopifnot_(
    !is.null(data[[outcome]]),
    "Can't find outcome variable {.var {outcome}} in {.arg data}."
  )
  error <- try_(match.arg(error, c("iid", "ar1")))
  stopifnot_(
    !inherits(error, "try-error"),
    "Argument {.arg error} must be either {.val iid} or {.val ar1}."
  )
  data <- data |>
    dplyr::arrange(.data[[unit]], .data[[time]]) |>
    droplevels()
  stopifnot_(
    all(table(data[[unit]], data[[time]]) == 1),
    "Data is not balanced."
  )
  treatment_table <- table(data[[treatment]], data[[unit]])
  treated_idx <- which(treatment_table[2, ] > 0)
  treated <- colnames(treatment_table)[treated_idx]
  N <- length(treated)
  stopifnot_(
    N >= 1L,
    "No treated units found in the data."
  )
  for (i in seq_along(treated)) {
    stopifnot_(
      !is.unsorted(data[[treatment]][data[[unit]] == treated[i]]),
      "There should be no gaps in the treatment. Check the treatment 
      variable for unit {.val {treated[i]}}."
    )
  }
  donors <- names(which(treatment_table[2, ] == 0))
  T_pre <- unname(treatment_table[1, treated])
  T_total <- length(unique(data[[time]]))
  stopifnot_(
    error == "iid" || min(T_pre) > 1L,
    "When {.arg error} is {.val ar1}, all treated units must have at least 
    two pre-treatment time points."
  )
  Y <- data |>
    dplyr::filter(.data[[unit]] %in% .env$treated) |>
    dplyr::pull(.data[[outcome]]) |>
    matrix(nrow = T_total)
  Z <- data |>
    dplyr::filter(.data[[unit]] %in% .env$donors) |>
    dplyr::pull(.data[[outcome]]) |>
    matrix(nrow = T_total)
  stopifnot_(
    all(!is.na(Z)),
    "Missing values in donor units are not supported."
  )
  missing_idx <- do.call(
    rbind,
    lapply(seq_along(T_pre), function(i) {
      miss <- which(is.na(Y[seq_len(T_pre[i]), i]))
      if (length(miss) > 0) {
        cbind(unit = i, time = miss)
      } else {
        NULL
      }
    })
  )
  if (is.null(missing_idx)) {
    missing_idx <- matrix(
      nrow = 0,
      ncol = 2,
      dimnames = list(NULL, c("unit", "time"))
    )
  }
  Y[is.na(Y)] <- 0
  colnames(Y) <- treated
  colnames(Z) <- donors
  J <- ncol(Z)
  beta_names <- gamma_names <- NULL
  X_y <- X_z <- X <- NULL
  tv_idx <- numeric(0)
  if (has_x) {
    X <- stats::model.matrix(formula, data = data)
    if (has_icpt) {
      X <- X[, -1L, drop = FALSE]
    }
    beta_names <- colnames(X)
    stopifnot_(
      nrow(X) == nrow(data),
      "Missing covariate values are not supported."
    )
    K <- ncol(X)
    n_units <- J + N
    X <- simplify2array(
      lapply(seq_len(K), \(k) matrix(X[, k], n_units, T_total, TRUE))
    )
    donor_idx <- which(!colnames(treatment_table) %in% treated)
    X_z <- X[donor_idx, , , drop = FALSE]
    X_y <- X[treated_idx, , , drop = FALSE]
    if (has_w) {
      gamma_names <- colnames(
        stats::model.matrix(parsed_formula$w_formula, data = data)
      )
      gamma_names <- setdiff(gamma_names, "(Intercept)")
      tv_idx <- match(gamma_names, beta_names)
      stopifnot_(
        !anyNA(tv_idx),
        c(
          "Column mismatch between time-varying and full predictor matrices.",
          i = "This shouldn't be possible, please report a bug on Github."
        )
      )
    }
  }
  setup <- dplyr::lst(
    outcome,
    treatment,
    treated,
    donors,
    unit,
    time,
    times = unique(data[[time]]),
    T_pre,
    T_total,
    has_icpt = has_icpt,
    has_x = has_x,
    has_w = has_w,
    has_ar1 = error == "ar1",
    missing_idx = missing_idx,
    predictors,
    beta_names,
    gamma_names,
    tv_idx,
    spline_df,
    knot_spacing,
    spline_type,
    noncentered_xi,
    prior_only
  )
  stan_args <- list(...)
  stan_args$chains <- stan_args$chains %||% 4L
  if (is.null(stan_args$control$adapt_delta) && has_w) {
    stan_args$control$adapt_delta <- 0.95
  }
  if (is.null(stan_args$iter) && is.null(stan_args$warmup)) {
    stan_args$iter <- 5000L
    stan_args$warmup <- 2500L
  }
  stan_args$object <- stanmodels$bscm
  if (is.null(stan_args$include)) {
    stan_args$include <- FALSE
    exclude_pars <- c(
      if (!has_icpt) "alpha",
      if (!has_x) "beta",
      if (!has_w) c("gamma", "kappa"),
      if (error != "ar1") "rho"
    )
    exclude_extras <- c("eta", "a", "xi")
    if (is.null(stan_args$pars)) {
      stan_args$pars <- c(exclude_pars, exclude_extras)
    } else {
      stan_args$pars <- union(stan_args$pars, exclude_extras)
    }
    setup$excluded_pars <- stan_args$pars
  } else {
    if (!stan_args$include) {
      setup$excluded_pars <- stan_args$pars
    }
  }
  if (setup$has_w) {
    spline_def <- build_spline(
      T_total,
      T_pre,
      spline_df,
      knot_spacing,
      spline_type,
      noncentered_xi
    )
  } else {
    spline_def <- NULL
  }
  descriptives <- compute_descriptives(Y, Z, T_pre, X_y, X_z, beta_names)
  priors <- define_priors(priors, descriptives, setup, spline_def)
  stan_args$data <- create_standata(
    setup,
    priors,
    Y,
    Z,
    X_y,
    X_z,
    spline_def,
    prior_only
  )
  if (is.null(stan_args$init)) {
    stan_args$init <- replicate(
      stan_args$chains,
      create_inits(stan_args$data, descriptives, spline_def),
      simplify = FALSE
    )
  }
  start_time <- proc.time()
  fit <- do.call(rstan::sampling, stan_args)
  stan_args$data$sample_y_rep <- compute_predictions
  gq <- rstan::gqs(
    stanmodels$generated_quantities,
    data = stan_args$data,
    draws = as.matrix(fit)
  )

  out <- list(
    stanfit = fit,
    y_mean = as.matrix(gq, "y_mean"),
    y_rep = if (compute_predictions) as.matrix(gq, "y_rep") else "Not sampled",
    data = if (save_data) data else NULL,
    setup = setup,
    priors = priors
  )
  class(out) <- "bscmfit"

  run_diags <- mcmc_diagnostics &&
    ndraws(out) > 50L &&
    !identical(stan_args$algorithm, "Fixed_param")
  if (run_diags) {
    out$converge <- check_mcmc_diagnostics.bscmfit(out)
  }
  out$elapsed_time <- list(
    total = proc.time() - start_time,
    sampling = rstan::get_elapsed_time(fit)
  )
  out$formula <- formula
  out$call <- match.call()
  out
}
