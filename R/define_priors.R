default_priors <- function(d, setup, spline_def) {
  N <- length(setup$treated)
  J <- length(setup$donors)
  pr_sigma <- half_normal_pr(signif(d$sd_e, 2))
  pr_omega <- logistic_normal_pr(
    default_zeta_prior(N, J),
    tilt = if (setup$marginalize_sigma) list(R = 3, scale = 2)
  )
  pr_alpha <- pr_beta <- pr_kappa <- pr_nu <- pr_rho <- NULL
  if (setup$has_icpt) {
    pr_alpha <- normal_pr(
      signif(d$mean_y, 2),
      signif(2 * d$sd_e, 2)
    )
  }
  if (setup$has_x) {
    sd_yx <- d$sd_y / d$sd_x
    pr_beta <- normal_pr(0, signif(2 * sd_yx, 2))
    if (setup$has_w) {
      pr_kappa <- half_normal_pr(signif(0.5 * sd_yx[setup$tv_idx], 2))
      if (setup$spline_type == "rw2") {
        c_nu <- 2 * sqrt(2) / (max(setup$T_pre) - 1)
        pr_nu <- normal_pr(0, signif(c_nu * sd_yx[setup$tv_idx], 2))
      }
    }
  }
  if (setup$has_ar1) {
    pr_rho <- beta_pr(rep(2, N), rep(2, N))
  }
  priors <- list(
    intercept = pr_alpha,
    coef = pr_beta,
    tv_sd = pr_kappa,
    tv_slope = pr_nu,
    sigma = pr_sigma,
    ar1 = pr_rho,
    weights = pr_omega
  )
  priors[lengths(priors) > 0]
}

#' Names of the model parameters in Stan corresponding to the role names
#'
#' The `priors` list and the user-facing output use the role names, while the
#' Stan models and the posterior draws use the names of the model parameters.
#' `concentration` and `tilt` are defined inside the prior of `weights`.
#'
#' @noRd
stan_parameter_names <- c(
  intercept = "alpha",
  coef = "beta",
  tv_sd = "kappa",
  tv_slope = "nu",
  sigma = "sigma",
  ar1 = "rho",
  weights = "omega",
  concentration = "zeta",
  tilt = "phi"
)

#' Default prior of the concentrations of the prior of the donor weights
#'
#' @param N \[`integer(1)`]\cr Number of treated units.
#' @param J \[`integer(1)`]\cr Number of donors.
#' @noRd
default_zeta_prior <- function(N, J) {
  a <- signif(sqrt(J), 2)
  b <- signif(0.5 * sqrt(J) / log(1 + J / 5), 2)
  inv_gamma_pr(rep(a, N), rep(b, N))
}

define_priors <- function(priors, descriptives, setup, spline_def) {
  defaults <- default_priors(descriptives, setup, spline_def)
  if (is.null(priors)) {
    return(defaults)
  }
  simplex_priors <- c(
    "dirichlet", "logistic_laplace", "logistic_normal"
  )
  bounded_priors <- c(
    "exponential", "gamma", "half_normal", "inv_gamma", "scaled_inv_chi",
    "sqrt_inv_gamma"
  )
  unbounded_priors <- c(
    "cauchy", "laplace", "normal", "student_t"
  )
  # unbounded can be used for bounded but not vice versa
  valid_priors <- list(
    intercept = unbounded_priors,
    coef = unbounded_priors,
    tv_sd = c(bounded_priors, unbounded_priors),
    tv_slope = unbounded_priors,
    sigma = c(bounded_priors, unbounded_priors),
    ar1 = "beta",
    weights = simplex_priors,
    concentration = c(bounded_priors, unbounded_priors)
  )

  invalid_priors <- setdiff(names(priors), names(defaults))
  stopifnot_(
    length(invalid_priors) == 0,
    c(
      "Invalid parameter {cli::qty(invalid_priors)}name{?s} in {.arg priors}:
      {.val {invalid_priors}}.",
      i = "Parameters of this model are {.val {names(defaults)}}."
    )
  )
  N <- length(setup$treated)
  for (nm in names(priors)) {
    pr <- priors[[nm]]
    dist <- pr$distribution
    stopifnot_(
      dist %in% valid_priors[[nm]],
      c(
        "Invalid prior family '{dist}' for parameter '{nm}'.",
        i = "Allowed families: {valid_priors[[nm]]}"
      )
    )
    if (nm == "weights" && !is.numeric(pr$concentration)) {
      zeta_prior <- pr$concentration %||% defaults$weights$concentration
      dist <- zeta_prior$distribution
      stopifnot_(
        dist %in% valid_priors$concentration,
        c(
          "Invalid prior family '{dist}' for the concentration of the prior
          of 'weights'.",
          i = "Allowed families: {valid_priors$concentration}"
        )
      )
      pr$concentration <- recycle_prior(zeta_prior, "concentration", N)
      pr$length <- N
    }
    priors[[nm]] <- recycle_prior(pr, nm, defaults[[nm]]$length %||% N)
  }
  defaults[names(priors)] <- priors
  defaults
}

#' Recycle a scalar prior to the required length
#'
#' @param pr \[`bscm_prior`]\cr Prior to recycle
#' @param nm \[`character(1)`]\cr Name of the parameter
#' @param n_ref \[`integer(1)`]\cr Required length of the prior.
#' @param call The calling environment reported in the error message.
#' @noRd
recycle_prior <- function(pr, nm, n_ref, call = rlang::caller_env()) {
  n_pr <- pr$length
  if (n_pr == n_ref) {
    return(pr)
  }
  stopifnot_(
    n_pr == 1,
    c(
      "Length mismatch for parameter '{nm}'.",
      i = "Got {n_pr}, expected {n_ref}.",
      i = "Only scalar priors (length = 1) can be recycled."
    ),
    call = call
  )
  pr$location <- rep(pr$location, n_ref)
  pr$scale <- rep(pr$scale, n_ref)
  pr$shape <- rep(pr$shape, n_ref)
  pr$rate <- rep(pr$rate, n_ref)
  pr$shape1 <- rep(pr$shape1, n_ref)
  pr$shape2 <- rep(pr$shape2, n_ref)
  pr$concentration <- rep(pr$concentration, n_ref)
  pr$df <- rep(pr$df, n_ref)
  pr$length <- n_ref
  pr
}
