#' Estimate Bayesian synthetic control model
#' 
#' Function `bscm` estimates a Bayesian synthetic control models of varying
#' types.
#' 
#' @details
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
#' The default priors for intercept \eqn{\alpha_i} and standard deviation 
#' \eqn{\sigma_i} of the noise term are 
#' \eqn{\alpha_i \sim N(m_i, s_i^2)} and
#' \eqn{\sigma_i \sim exponential(l_i)},
#' where by default `m_i = mean(y_i)`, `s_i = 2 * sd(y_i)`,
#' `l_i = 1 / sd(y_i)`, with `y_i` being the vector of pretreatment
#' outcomes of treated unit \eqn{i}.
#' 
#' The default prior for the coefficient \eqn{\beta_k} is  
#' \eqn{\beta_k \sim N(0, s_k^2)}, where \eqn{s_k = 2 * s / sd(x_k)} and s
#' is the mean of the variances of the treated and donors' outcomes in the 
#' pretreatment period, and \eqn{x_k} is the vector of values of the k-th 
#' covariate in the pretreatment period over all.
#' 
#' The prior for the weight vector \eqn{\omega} is defined as a symmetric 
#' Dirichlet distribution with a concentration parameter \eqn{\kappa}. Value 
#' \eqn{\kappa = 1} corresponds to uniform distribution on a
#' simplex. When \eqn{\kappa < 1}, the prior prefers sparse weight vectors, 
#' while values larger than one put more prior probability mass on solutions 
#' where \eqn{\omega} is dense. The default is \eqn{\kappa = 0.5}.
#' 
#' Alternatively, argument `effective_donors` can be used to define 
#' \eqn{\kappa} using the approximation of expected number of 
#' effective donors, defined as \eqn{1 / \sum_{j=1}^J(\omega_j^2)}. More 
#' specifically, `kappa = (effective_donors - 1) / (J - effective_donors)` 
#' where \eqn{J} is the total number of donors. Default \eqn{\kappa = 0.5} 
#' therefore corresponds to prior where the effective number donors is
#' \eqn{(J + 2) / 3}.
#' 
#' The output of `bscm()` contains posterior samples of various derived 
#' quantities such as effect estimates. To access these after model
#' estimation, use methods such as [treatment_effect()] [coef()], and 
#' [summary.bscmfit()] on the output of `bscm()`. 
#' 
#' @param formula \[`formula`]\cr The model formula containing the outcome 
#' variable on the left-hand side and optional time-varying predictors on
#' the right-hand side (RHS) of `~`. In case of no predictors, use
#' `outcome ~ 1` or `outcome ~ 0`. In the former case, as well by the
#' default when using predictors (e.g., `outcome ~ x + z`), the model
#' includes intercept term. Intercept can be omitted by using `0` in the
#' RHS, e.g., `outcome ~ 0` or `outcome ~ 0 + x + z` (equivalently, you
#' can use `-1` in place of `0`). Formula does not need to contain the
#' variable defining the treatment, which is defined separately using the
#' argument `treatment`. In case the variable is present also in the
#' formula, it is automatically removed.
#' 
#' To specify predictors with time-varying coefficients, wrap them in
#' `tv()`, e.g., `outcome ~ x + tv(~ z)` defines a model where `x` has
#' a time-constant coefficient and `z` has a time-varying coefficient.
#' Terms inside `tv()` are automatically also included in the
#' time-constant part of the model, since the time-varying coefficients
#' are defined to have zero mean. Possible time-varying intercept is omitted as 
#' it would cancel out in the linear predictor.
#' 
#' @param data \[`data.frame` or an object coercible to one]\cr
#'   The long format data that contains the model variables.
#' @param treatment  \[`character(1)`]\cr Name of the treatment indicator 
#'   variable in `data`.
#' @param time \[`character(1)`]\cr Name of the time index variable in
#'   `data`.
#' @param unit \[`character(1)`]\cr Name of the variable in `data`
#'   identifying different units.
#' @param priors \[`list()` or `character(1)`]\cr List of prior definitions 
#'   or `"default"` which is a default and only supported option at the
#'   moment. See details on the prior definitions.
#' @param kappa \[`numeric(1)`]\cr A positive number defining the
#'   concentration parameter \eqn{\kappa} of symmetric Dirichlet prior of
#'   the donor weights. Ignored if the argument `effective_donors` is
#'   non-NULL. Note that small values of \eqn{\kappa}  can lead to
#'   divergences in sampling.
#' @param effective_donors \[`integer(1)`]\cr Integer for alternative 
#'   definition of the prior for the donor weights. See details.
#' @param mcmc_diagnostics \[`logical(1)`]\cr If `TRUE` (the default), the 
#'   output of [bscm()] includes the results of MCMC diagnostics checks
#'   performed by [check_mcmc_diagnostics.bscmfit()]. Note that regardless
#'   of the value of `mcmc_diagnostics`, [rstan::sampling()] can still
#'   generate warnings regarding convergence and other sampling issues.
#' @param save_data \[`logical(1)`]\cr If `TRUE` (the default), `bscmfit` 
#'   object returned by  [bscm()] includes the input `data.frame` 
#'   (argument `data`), after dropping unused factor levels and potentially 
#'   rearranging data by `unit` and `time` variables.
#' @param ... Additional parameters passed on to [rstan::sampling()] to
#'   adjust the sampling options, for example `iter` and `chains`. Note that
#'   defaults `iter = 5000` and `warmup = 2500` differ from the defaults of 
#'   [rstan::sampling()] (which are 2000 and 1000 respectively).   
#' @return An object of class `bscmfit`.
#' @export
#' @seealso [summary.bscmfit()], [as_draws.bscmfit()], [rstan::sampling()].
#' @examples
#' fit <- bscm(
#'   y ~ 1, single_treated, "treatment", "time", "id", 
#'   chains = 1, cores = 1, refresh = 0
#' )
#' fit
#' 
bscm <- function(formula, data, treatment, time = "time", unit = "id",
                 priors = "default", kappa = 0.5, effective_donors = NULL,
                 mcmc_diagnostics = TRUE, save_data = TRUE, ...) {
  
  # local function for creating input data to Stan
  create_standata <- \() {
    s_y <- pmax(1, sd_y)
    standata <- list(
      T = T_total, T_pre = array(T_pre), N = N, J = J, 
      y = t(Y), Z = Z,
      pr_rate_sigma = array(1 / s_y),
      pr_mean_intercept = array(mean_y), 
      pr_sd_intercept = array(2 * s_y),
      kappa = kappa
    )
    if (has_x) {
      s_x <- pmax(1, sd_x)
      s_yz <- max(1, sqrt(stats::median(c(sd_y, sd_z)^2)))
      pr_sd_beta <- 2 * s_yz / s_x
      standata <- c(
        standata, 
        list(
          K = K, X_z = aperm(X_z, c(3, 2, 1)), X_y = X_y,
          pr_mean_beta = array(0, K), pr_sd_beta = array(pr_sd_beta)
        )
      )
      if (has_w) {
        s_w <- pmax(1, sd_w)
        pr_rate_tau <- 5 * s_w / s_yz
        standata <- c(
          standata, 
          list(
            L = L, tv_idx = array(tv_idx),
            pr_rate_tau = array(pr_rate_tau)
          )
        )
      }
    }
    standata
  }
  
  # local function to create initial values for Stan
  create_inits <- \() {
    s_y <- pmax(1, sd_y)
    inits <- list(
      omega = matrix(1 / J, N, J),
      sigma = array(stats::runif(N, 0.5 * s_y, 2 * s_y))
    )
    if (has_icpt) {
      inits$a <- array(stats::rnorm(N, mean_y, 0.5 * s_y))
    }
    if (has_x) {
      s_x <- pmax(1, sd_x)
      s_yz <- max(1, sqrt(mean(c(sd_y, sd_z)^2)))
      pr_sd_beta <- 2 * s_yz / s_x
      inits$beta <- array(stats::rnorm(K, 0, 0.1 * pr_sd_beta))
      if (has_w) {
        inits$gamma_raw <- matrix(0, T_total, L)
        inits$tau <- array(stats::runif(L, 0.05, 0.1))
      }
    }
    inits
  }
  
  # start the actual function body
  tol <- sqrt(.Machine$double.eps)
  check_bscm_arguments(
    formula, data, treatment, time, unit, 
    priors, kappa, effective_donors, save_data
  )
  outcome <- get_outcome(formula)
  parsed_formula <- parse_bscm_formula(formula)
  formula <- parsed_formula$x_formula
  has_icpt <- parsed_formula$icpt
  predictors <- parsed_formula$predictors
  has_x <- length(predictors > 0)
  has_w <- length(parsed_formula$w_terms > 0)
  stopifnot_(
    !is.null(data[[outcome]]),
    "Can't find outcome variable {.var {outcome}} in {.arg data}."
  )
  data <- data |>
    arrange(.data[[unit]], .data[[time]]) |> 
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
  Y <- data |> 
    filter(.data[[unit]] %in% .env$treated) |> 
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
  mean_y <- vapply(
    seq_len(N), \(i) mean(Y[seq_len(T_pre[i]), i]), numeric(1)
  )
  sd_y <- vapply(
    seq_len(N), \(i) sd(Y[seq_len(T_pre[i]), i]), numeric(1)
  )
  stopifnot_(
    all(sd_y > tol),
    "Outcome variable cannot be constant in the pre-treatment period. 
      Found `sd(y) < sqrt(.Machine$double.eps)`."
  )
  min_T_pre <- min(T_pre)
  sd_z <- apply(Z[seq_len(min_T_pre), , drop = FALSE], 2, sd)
  mean_z <- colMeans(Z[seq_len(min_T_pre), , drop = FALSE])
  stopifnot_(
    any(sd_z > tol),
    "Outcome variable cannot be constant in the pre-treatment period. 
      Found `sd(z) < sqrt(.Machine$double.eps)`."
  )
  beta_names <- gamma_names <- NULL
  if (has_x) {
    X <- stats::model.matrix(formula, data = data)
    if (has_icpt) X <- X[, -1L, drop = FALSE]
    beta_names <- colnames(X)
    stopifnot_(
      nrow(X) == nrow(data),
      "Missing values are not supported."
    )
    K <- ncol(X)
    n_units <- J + N
    X <- simplify2array(
      lapply(seq_len(K), \(k) matrix(X[, k], n_units, T_total, TRUE))
    )
    sd_x_by_unit <- apply(
      X[, seq_len(min_T_pre), , drop = FALSE], c(1, 3), sd
    )
    constant_sd <- which(
      stats::setNames(apply(sd_x_by_unit, 2, max) < tol, beta_names)
    )
    warnifnot_(
      length(constant_sd) == 0 || !has_icpt,
      c(
        "Model has unit-specific intercepts and predictors which do not vary in 
        the pre-treatment period for any units.",
        i = "Found {?a/} constant predictor{?s} {names(constant_sd)}."
      )
    )
    sd_x <- apply(sd_x_by_unit, 2, stats::median)
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
        "Column mismatch between time-varying and full predictor matrices."
      )
      L <- length(tv_idx)
      sd_w <- sd_x[tv_idx]
    }
  }
  if (!is.null(effective_donors)) {
    stopifnot_(
      effective_donors >= 2 && effective_donors < J - 1,
      "Argument {.arg effective_donors} should be between 2 and {J - 1}."
    )
    kappa <- (effective_donors - 1) / (J - effective_donors)
  }
  icpt <- ifelse(has_icpt, "int", "noint")
  x <- ifelse(has_x, "x", "nox")
  effect <- if (has_w) {
    "varying"
  } else if (has_x) {
    "const"
  } else {
    "none"
  }
  model_type <- paste("bscm", icpt, x, effect, sep = "_")
  
  stan_args <- list(...)
  stan_args$chains  <- stan_args$chains %||% 4L
  if (is.null(stan_args$iter) && is.null(stan_args$warmup)) {
    stan_args$iter <- 5000L
    stan_args$warmup <- 2500L
  }
  stan_args$data <- create_standata()
  if (is.null(stan_args$init)) {
    stan_args$init <- replicate(
      stan_args$chains, 
      create_inits(),
      simplify = FALSE
    )
  }
  stan_args$object <- stanmodels[[model_type]]
  stan_args$pars <- c(
    "omega_raw", 
    if (has_icpt) "a",
    if (has_w) "gamma_raw"
  )
  stan_args$include <- FALSE
  start_time <- proc.time()
  fit <- do.call(sampling, stan_args)
  out <- list(stanfit = fit)
  if (save_data) out$data <- data
  times <- unique(data[[time]])
  priors <- stan_args$data[substr(names(stan_args$data), 1, 3) %in% c("pr_")]
  out$setup <- dplyr::lst(
    outcome, treatment, treated, donors, unit, time, times, 
    T_pre, T_total, 
    has_intercept = has_icpt,
    predictors, beta_names, gamma_names, kappa,
    model_type, priors
  )
  class(out) <- "bscmfit"
  if (mcmc_diagnostics && !identical(stan_args$algorithm, "Fixed_param")) {
    out$converge <- check_mcmc_diagnostics.bscmfit(out)
  }
  out$elapsed_time <- list(
    total = proc.time() - start_time,
    sampling = rstan::get_elapsed_time(fit)
  )
  out$call <- match.call()
  out
}
