#' Estimate Bayesian synthetic control model
#'
#' Function `bscm` estimates a Bayesian synthetic control models of varying
#' types.
#'
#' @details
#'
#' The prior for the weight vector \eqn{\omega} is controlled by the
#' `omega_prior` argument. Two families are supported:
#'
#' - Logistic normal ([logistic_normal()]): \eqn{\omega =
#'   \textrm{softmax}(\eta)} with \eqn{\eta \sim N(0, \kappa^2 I)} constrained
#'   to sum to zero. Larger \eqn{\kappa} induces sparser weights.
#' - Symmetric Dirichlet ([dirichlet()]): \eqn{\omega \sim
#'   \textrm{Dirichlet}(\kappa, \ldots, \kappa)}. Values \eqn{\kappa < 1}
#'   concentrate weight on few donors while \eqn{\kappa > 1} pulls encourages
#'   more uniform weights. The default prior is
#'   \eqn{Dirichlet(\kappa = 1)} which corresponds to uniform prior over
#'   probability simplices.
#'
#' You should test different values to assess sensitivity of results,
#' and potentially run [bscm::loo()] or [bscm::lfo()] for cross-validation
#' based selection (although it can be inconclusive for small and moderate
#' number of pre-treatment time points and/or donors).
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
#' \eqn{\beta_k \sim N(0, s_k^2)}, where \eqn{s_k = 2s_y / s_{x,k}} and
#' \eqn{s_y} and  \eqn{s_{x,k}} are the pre-treatment median standard
#' deviations of the outcomes for treated units and covariates \eqn{x_k} for
#' all units.
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
#'   moment for parameters other than weight vector \eqn{\omega}. See details.
#' @param omega_prior \[`omega_prior`]\cr Prior for the donor weight vector
#'   \eqn{\omega}, created by [logistic_normal()] or [dirichlet()], where both
#'   constructors take argument `kappa` which defines the scale and
#'   concentration parameter of the corresponding distribution.
#'   Defaults to `dirichlet(kappa = 1)`. See details.
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
#' @param no_donors \[`logical(1)`]\cr Should donors be ignored? Default is
#'   `FALSE`, but if set to `TRUE`, instead of BSCM, a simple Bayesian
#'   linear regression model based on `formula` is estimated. This is mainly
#'   for the purposes of assessing the relative the performance of BSCM in
#'   case where convex hull assumption of SCM does not hold.
#' @param ... Additional parameters passed on to [rstan::sampling()] to
#'   adjust the sampling options, for example `iter` and `chains`. Note that
#'   defaults `iter = 5000` and `warmup = 2500` differ from the defaults of
#'   [rstan::sampling()] (which are 2000 and 1000 respectively).
#' @return An object of class `bscmfit`.
#' @export
#' @seealso [summary.bscmfit()], [as_draws.bscmfit()], [rstan::sampling()].
#' @examples
#' # skip diagnostics and use small number of iterations for CRAN checks
#' fit <- bscm(
#'   y ~ 1, single_treated, "treatment", "time", "id",
#'   omega_prior = dirichlet(0.5),
#'   chains = 1, cores = 1, refresh = 0, iter = 1000,
#'   mcmc_diagnostics = FALSE
#' )
#' fit
bscm <- function(
  formula,
  data,
  treatment,
  time = "time",
  unit = "id",
  omega_prior = dirichlet(kappa = 1),
  mcmc_diagnostics = TRUE,
  save_data = TRUE,
  priors = "default",
  compute_predictions = TRUE,
  no_donors = FALSE,
  ...
) {
  check_bscm_arguments(
    formula,
    data,
    treatment,
    time,
    unit,
    omega_prior,
    mcmc_diagnostics,
    save_data,
    priors,
    compute_predictions,
    no_donors
  )
  outcome <- get_outcome(formula)
  parsed_formula <- parse_bscm_formula(formula)
  formula <- parsed_formula$x_formula
  has_icpt <- parsed_formula$icpt
  predictors <- parsed_formula$predictors
  has_x <- length(predictors > 0)
  has_w <- length(parsed_formula$w_terms > 0)
  stopifnot_(
    !no_donors || (no_donors && has_icpt),
    "Argument {.arg formula} should include intercept 
    when {.arg no_donors} is `TRUE`."
  )
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
  beta_names <- gamma_names <- NULL
  if (has_x) {
    X <- stats::model.matrix(formula, data = data)
    if (has_icpt) {
      X <- X[, -1L, drop = FALSE]
    }
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
      X[, seq_len(min(T_pre)), , drop = FALSE],
      c(1, 3),
      sd
    )
    constant_sd <- which(
      stats::setNames(
        apply(sd_x_by_unit, 2, max) < sqrt(.Machine$double.eps),
        beta_names
      )
    )
    warnifnot_(
      length(constant_sd) == 0 || !has_icpt,
      c(
        "Model has unit-specific intercepts and predictors which do not vary in 
        the pre-treatment period for any units.",
        i = "Found {?a/} constant predictor{?s} {names(constant_sd)}."
      )
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
        "Column mismatch between time-varying and full predictor matrices."
      )
      L <- length(tv_idx)
    }
  }

  icpt <- ifelse(has_icpt, "a1", "a0")
  x <- ifelse(has_x, "x1", "x0")
  w <- ifelse(has_w, "w1", "w0")
  omega_prior_type <- omega_prior$distribution
  o <- ifelse(omega_prior_type == "dirichlet", "dr", "ln")
  omega <- ifelse(no_donors, "no", o)
  model_type <- paste(
    c("bscm", icpt, x, w, omega),
    collapse = "_"
  )

  stan_args <- list(...)
  stan_args$chains <- stan_args$chains %||% 4L
  if (is.null(stan_args$control$adapt_delta)) {
    stan_args$control$adapt_delta <- 0.95
  }
  if (is.null(stan_args$iter) && is.null(stan_args$warmup)) {
    stan_args$iter <- 5000L
    stan_args$warmup <- 2500L
  }
  stan_args$object <- stanmodels[[model_type]]
  if (is.null(stan_args$pars) && is.null(stan_args$include)) {
    exclude_these <- c(
      if (omega_prior_type == "logistic_normal") "eta",
      if (has_icpt) "a",
      if (has_w) "gamma_raw"
    )
    if (length(exclude_these) > 0L) {
      stan_args$pars <- exclude_these
      stan_args$include <- FALSE
    }
  }

  input_stats <- bscm_stats(Y, Z, T_pre, X = if (has_x) X)
  stan_args$data <- create_standata(
    input_stats,
    T_pre,
    Y,
    Z,
    has_icpt,
    omega_prior$kappa,
    X_y = if (has_x) X_y,
    X_z = if (has_x) X_z,
    tv_idx = if (has_w) tv_idx,
    cv = as.integer(!compute_predictions)
  )
  if (is.null(stan_args$init)) {
    stan_args$init <- replicate(
      stan_args$chains,
      create_inits(stan_args$data, omega_prior),
      simplify = FALSE
    )
  }

  start_time <- proc.time()
  fit <- do.call(sampling, stan_args)
  out <- list(stanfit = fit)
  if (save_data) {
    out$data <- data
  }
  times <- unique(data[[time]])
  priors <- stan_args$data[startsWith(names(stan_args$data), "pr_")]
  out$setup <- dplyr::lst(
    outcome,
    treatment,
    treated,
    donors,
    unit,
    time,
    times,
    T_pre,
    T_total,
    has_intercept = has_icpt,
    predictors,
    beta_names,
    gamma_names,
    model_type,
    priors,
    omega_prior
  )
  class(out) <- "bscmfit"
  if (mcmc_diagnostics && !identical(stan_args$algorithm, "Fixed_param")) {
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
