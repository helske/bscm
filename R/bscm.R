#' Estimate Bayesian synthetic control model
#' 
#' Function `bscm` estimates a Bayesian synthetic control models of varying types.
#' 
#' @details
#' 
#' When model contains covariates \eqn{X}, their effect is subtracted from 
#' donors, i.e., for treated unit 
#' \eqn{y \sim N(\alpha + X\beta + Z^\ast\omega, \sigma^2)}, 
#' where column \eqn{j} of \eqn{Z^\ast} is \eqn{z_j - (\alpha_j + X_j \beta)}, 
#' i.e., the original donor outcome vector \eqn{z_j} minus the donor-specific 
#' intercept \eqn{\alpha_j} and the effect of covariates \eqn{X_j \beta} for 
#' that donor. Note that \eqn{\beta} is common across donors and the treated.
#' 
#' The default priors for intercept \eqn{\alpha} and standard deviation 
#' \eqn{\sigma} of the noise term are 
#' \eqn{\alpha \sim N(m, s^2)} and \eqn{\sigma \sim exponential(l)},
#' where by default `m = mean(y)`, `s = 2 * sd(y)`, `l = 1 / sd(y)`,
#' with `y` being the vector of pretreatment outcomes. Same priors are used for 
#' the donors as well, with `m` and `s` defined separately for each donor using 
#' the pretreatment outcomes of that donor.
#' 
#' The default prior for the coefficient \eqn{\beta_k} is  
#' \eqn{\beta_k \sim N(0, s_k^2)}, where \eqn{s_k = 2 * s / sd(x_k)} and s is 
#' the mean of the variances of the treated and donors' outcomes in the 
#' pretreatment period, and \eqn{x_k} is the vector of values of the k-th 
#' covariate in the pretreatment period over all.
#' 
#' The prior for the weight vector \eqn{\omega} is defined as a symmetric 
#' Dirichlet distribution with a concentration parameter \eqn{\kappa}. By 
#' default, \eqn{\kappa = 1}, corresponding to uniform distribution on a simplex. 
#' Alternatively, argument `effective_donors` can be used to define 
#' \eqn{\kappa} using the approximation of expected number of 
#' effective donors, defined as \eqn{1 / \sum(w^2)}. More specifically, 
#' `kappa = (effective_donors - 1) / (J - effective_donors)` where \eqn{J} is 
#' the total number of donors. Default \eqn{\kappa = 1} therefore corresponds 
#' to prior where the effective number donors is \eqn{(J + 1) / 2}.
#' 
#' The output of `bscm()` contains posterior samples of various derived 
#' quantities such as effect estimates. To access these after model estimation,
#' use [summary.bscmfit()] method on the output of `bscm(). 
#' 
#' @param formula \[`formula`]\cr The model formula containing the outcome 
#' variable on the left-hand side and optional time-varying predictors on the 
#' right-hand side (RHS) of `~`. In case of no predictors, use `outcome ~ 1` or 
#' `outcome ~ 0`. In the former case, as well by the default when using 
#' predictors (e.g., `outcome ~ x + z`), the model includes intercept term. 
#' Intercept can be omitted by using `0` in the RHS, e.g.,
#' `outcome ~ 0` or `outcome ~ 0 + x + z` (equivalently, you can use `-1` in 
#' place of `0`). Formula does not need to contain the variable defining the 
#' treatment, which is defined separately using the argument `treatment`. In 
#' case the variable is present also in the formula, it is automatically 
#' removed.
#' 
#' @param data \[`data.frame` or an object coercible to one]\cr
#'   The long format data that contains the model variables.
#' @param treatment  \[`character(1)`]\cr Name of the treatment indicator 
#' variable in `data`.
#' @param time \[`character(1)`]\cr Name of the time index variable in `data`.
#' @param unit \[`character(1)`]\cr Name of the variable in `data` identifying 
#' different units.
#' @param time_varying_effects \[`logical(1)`]\cr Should the
#' regression coefficients vary over time? The default is `FALSE`, while 
#' `TRUE` defines the regression coefficients as penalized splines. See details.
#' @param priors \[`list()` or `character(1)`]\cr List of prior definitions 
#' or `"default"` which is a default and only supported option at the moment. 
#' See details on the prior definitions.
#' @param kappa \[`numeric(1)`]\cr A positive number defining the concentration 
#' parameter of symmetric Dirichlet prior of the donor weights. Ignored if the 
#' argument `effective_donors` is non-NULL.
#' @param effective_donors \[`integer(1)`]\cr Integer for alternative 
#' definition of the prior for the donor weights. See details.
#' @param mcmc_diagnostics \[`logical(1)`]\cr If `TRUE` (the default), the 
#' output of [bscm()] includes the results of MCMC diagnostics checks performed 
#' by [check_mcmc_diagnostics.bscmfit()]. Note that regardless of 
#' the value of `mcmc_diagnostics`, [rstan::sampling()] can still generate 
#' warnings regarding convergence and other sampling issues.
#' @param save_data \[`logical(1)`]\cr If `TRUE` (the default), `bscmfit` 
#' object returned by  [bscm()] includes the input `data.frame` 
#' (argument `data`), after dropping unused factor levels and potentially 
#' rearranging data by `unit` and `time` variables.
#' @param ... Additional parameters passed on to [rstan::sampling()] to adjust 
#' the sampling options.
#' @return An object of class `bscmfit`.
#' @export
#' @seealso [summary.bscmfit()], [as_draws.bscmfit()], [rstan::sampling()].
#' @examples
#' fit <- bscm(
#'   y ~ 1, simulated_data, "treatment", "time", "id", 
#'   chains = 1, cores = 1, refresh = 0
#' )
#' fit
#' 
bscm <- function(formula, data, treatment, time = "time", unit = "id",
                 time_varying_effects = FALSE,
                 priors = "default", kappa = 1, effective_donors = NULL,
                 mcmc_diagnostics = TRUE,
                 save_data = TRUE, ...) {
  
  # local function for creating input data to Stan
  create_standata <- \() {
    s_y <- max(1, sd_y)
    standata <- list(
      T = T_total, T_pre = T_pre, J = J, y = c(Y), Z = Z, kappa = kappa,
      pr_rate_sigma = 1 / s_y,
      pr_mean_intercept = mean_y, pr_sd_intercept = 2 * s_y
    )
    if (has_predictors) {
      s_x <- pmax(1, sd_x)
      s_yz <- max(1, sqrt(mean(c(sd_y, sd_z)^2)))
      pr_sd_coef <- 2 * s_yz / s_x
      standata <- c(
        standata, 
        list(
          K = K, X_z = X_z, X_y = X_y,
          pr_mean_coef = array(0, K), pr_sd_coef = array(pr_sd_coef, K),
          inv_sd_x = array(1 / sd_x, K)
        )
      )
      if (effect == "varying") {
        standata$pr_rate_tau <- 1
      }
    }
    standata
  }
  
  # local function to create initial values for Stan
  create_inits <- \() {
    s_y <- max(1, sd_y)
    inits <- list(
      omega_raw = rep(1, J),
      sigma = stats::runif(1, 0.5 * s_y, 2 * s_y)
    )
    if (has_intercept) {
      inits$a <- stats::rnorm(1, mean_y, 0.1)
    }
    if (has_predictors) {
      s_x <- pmax(1, sd_x)
      s_yz <- max(1, sqrt(mean(c(sd_y, sd_z)^2)))
      pr_sd_coef <- 2 * s_yz / s_x
      inits$beta <- array(stats::rnorm(K, 0, 0.1 * pr_sd_coef), K)
      if (time_varying_effects) {
        inits$tau <- stats::runif(K, 0.5, 2)
      }
    }
    inits
  }
  
  # start the actual function body
  tol <- sqrt(.Machine$double.eps)
  check_bscm_arguments(
    formula, data, treatment, time, unit, time_varying_effects, 
    priors, kappa, effective_donors, save_data
  )
  outcome <- get_outcome(formula)
  stopifnot_(
    !is.null(data[[outcome]]),
    "Can't find outcome variable {.var {outcome}} in {.arg data}."
  )
  data <- data |>
    arrange(.data[[unit]], .data[[time]]) |> 
    droplevels()
  stopifnot_(
    all(table(data[[unit]], data[[time]]) == 1),
    "Data is not balanced. All units should have equal number of time points."
  )
  treatment_table <- table(data[[treatment]], data[[unit]])
  treated_idx <- which(treatment_table[2, ] > 0)
  treated <- colnames(treatment_table)[treated_idx]
  N <- length(treated)
  stopifnot_(
    identical(N, 1L),
    "Only the case of a single treated unit is currently implemented."
  )
  stopifnot_(
    !is.unsorted(data[[treatment]][data[[unit]] == treated]),
    "There should be no gaps in the treatment. Check the treatment variable."
  )
  donors <- names(which(treatment_table[2, ] == 0))
  T_pre <- treatment_table[1, treated]
  T_total <- length(unique(data[[time]]))
  Y <- data |> 
    filter(.data[[unit]] == .env$treated) |> 
    pull(.data[[outcome]]) |> 
    matrix(nrow = T_total)
  Z <- data |> 
    filter(.data[[unit]] %in% .env$donors) |> 
    pull(.data[[outcome]]) |> 
    matrix(nrow = T_total)
  stopifnot_(
    all(!is.na(Y)) && all(!is.na(Z)),
    "Missing values are not supported."
  )
  J <- ncol(Z)
  mean_y <- mean(Y[1:T_pre])
  sd_y <- sd(Y[1:T_pre])
  stopifnot_(
    sd_y > tol,
    "Outcome variable cannot be constant in the pre-treatment period. 
      Found `sd(y) < sqrt(.Machine$double.eps)`."
  )
  sd_z <- apply(Z[1:T_pre, , drop = FALSE], 2, sd)
  mean_z <- colMeans(Z[1:T_pre, , drop = FALSE])
  stopifnot_(
    any(sd_z > tol),
    "Outcome variable cannot be constant in the pre-treatment period. 
      Found `sd(z) < sqrt(.Machine$double.eps)`."
  )
  has_intercept <- as.logical(attr(stats::terms(formula), "intercept"))
  predictors <- all.vars(formula[-2L])
  has_predictors <- length(predictors) > 0L
  if (!has_predictors) time_varying_effects <- FALSE
  coef_names <- NULL
  if (has_predictors) {
    X <- stats::model.matrix(formula, data = data)
    if (has_intercept) X <- X[, -1L, drop = FALSE]
    coef_names <- colnames(X)
    stopifnot_(
      nrow(X) == nrow(data),
      "Missing values are not supported."
    )
    K <- ncol(X)
    X <- simplify2array(
      lapply(seq_len(K), \(k) matrix(X[, k], J + 1, T_total, TRUE))
    )
    sd_x_by_unit <- apply(X[, 1:T_pre, , drop = FALSE], c(1, 3), sd)
    constant_sd <- which(
      stats::setNames(apply(sd_x_by_unit, 2, max) < tol, coef_names)
    )
    stopifnot_(
      length(constant_sd) == 0,
      c(
        "Predictor variables cannot be constant in the pre-treatment period for 
      all units.",
        i = "Found {?a/} constant predictor(?s) {names(constant_sd)}."
      )
    )
    sd_x <- apply(sd_x_by_unit, 2, stats::median)
    X_z <- X[-treated_idx, , , drop = FALSE]
    X_y <- X[treated_idx, , , drop = FALSE]
    dim(X_y) <- c(T_total, K)
  }
  if (!is.null(effective_donors)) {
    stopifnot_(
      effective_donors >= 2 && effective_donors < J - 1,
      "Argument {.arg effective_donors} should be between 2 and {J - 1}."
    )
    kappa <- (effective_donors - 1) / (J - effective_donors)
  }
  warnifnot_(
    kappa > 0.1,
    "{.arg effective_donors} = {effective_donors} implies a Dirichlet prior 
        with concentration parameter {kappa}. Such a small value can lead 
        to divergences in sampling. If this happens, try increasing 
        {.arg control$adapt_delta} closer to 1` or increase 
        {.arg effective_donors}."
  )
  icpt <- ifelse(has_intercept, "int", "noint")
  x <- ifelse(has_predictors, "x", "nox")
  effect <- case_when(
    !has_predictors ~ "none",
    time_varying_effects ~ "varying",
    .default = "const"
  )
  model_type <- paste("bscm", icpt, x, effect, sep = "_")
  
  stan_args <- list(...)
  stan_args$chains  <- stan_args$chains %||% 4L
  stan_args$data <- create_standata()
  if (is.null(stan_args$init)) {
    stan_args$init <- replicate(
      stan_args$chains, 
      create_inits(),
      simplify = FALSE
    )
  }
  stan_args$object <- stanmodels[[model_type]]
  stan_args$pars <- c("omega_raw", if (has_intercept) "a")
  stan_args$include <- FALSE
  start_time <- proc.time()
  fit <- do.call(sampling, stan_args)
  out <- list(stanfit = fit)
  if (save_data) out$data <- data
  times <- unique(data[[time]])
  priors <- stan_args$data[substr(names(stan_args$data), 1, 3) %in% c("pr_")]
  out$setup <- dplyr::lst(
    outcome, treatment, treated, donors, unit, time, times, 
    T_pre, T_total, has_intercept, predictors, coef_names, kappa,
    model_type, time_varying_effects, priors
  )
  class(out) <- "bscmfit"
  if (mcmc_diagnostics) {
    out$converge <- check_mcmc_diagnostics.bscmfit(out)
  }
  out$elapsed_time <- list(
    total = proc.time() - start_time,
    sampling = rstan::get_elapsed_time(fit)
  )
  out$call <- match.call()
  out
}
