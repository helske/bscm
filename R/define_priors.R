default_priors <- function(d, setup, spline_def) {
  N <- length(setup$treated)
  J <- length(setup$donors)
  pr_sigma <- half_normal_pr(signif(0.5 * d$sd_e, 2))
  pr_omega <- logistic_normal_pr(rep(NA, N))
  pr_zeta <- inv_gamma_pr(rep(sqrt(J), N), rep(0.5 * sqrt(J) / log(1 + J / 5), N))
  
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
        pr_nu <- normal_pr(0, signif(2 * sd_yx[setup$tv_idx] / d$sd_t, 2))
      }
    }
  }
  if (setup$has_ar1) {
    pr_rho <- beta_pr(rep(2, N), rep(2, N))
  }
  priors <- list(
    sigma = pr_sigma,
    omega = pr_omega,
    zeta = pr_zeta,
    alpha = pr_alpha,
    beta = pr_beta,
    kappa = pr_kappa,
    nu = pr_nu,
    rho = pr_rho
  )
  priors[lengths(priors) > 0]
}

#' Plain-text synonyms for the names of the `priors` list
#'
#' @noRd
prior_aliases <- c(
  intercept = "alpha",
  coef = "beta",
  scale_gamma = "kappa",
  slope_gamma = "nu",
  ar1 = "rho",
  weights = "omega",
  concentration = "zeta",
  sd = "sigma"
)

#' Rename the elements of `priors` given by a synonym to the canonical name
#'
#' @param priors \[`list()`]\cr The `priors` argument of [bscm()].
#' @param call The calling environment reported in the error message.
#' @noRd
resolve_prior_aliases <- function(priors, call = rlang::caller_env()) {
  idx <- match(names(priors), names(prior_aliases))
  if (all(is.na(idx))) {
    return(priors)
  }
  canonical <- names(priors)
  canonical[!is.na(idx)] <- prior_aliases[idx[!is.na(idx)]]
  duplicated_priors <- unique(canonical[duplicated(canonical)])
  stopifnot_(
    length(duplicated_priors) == 0L,
    c(
      "{cli::qty(duplicated_priors)}Parameter{?s} {.val {duplicated_priors}}
      {?is/are} given multiple priors in {.arg priors}.",
      i = "A parameter can be named either by its canonical name or by a
      synonym, but not by both."
    ),
    call = call
  )
  stats::setNames(priors, canonical)
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
    sigma = c(bounded_priors, unbounded_priors),
    omega = simplex_priors,
    zeta = c(bounded_priors, unbounded_priors),
    alpha = unbounded_priors,
    beta = unbounded_priors,
    kappa = c(bounded_priors, unbounded_priors),
    nu = unbounded_priors,
    rho = "beta"
  )
  omega_defined <- "omega" %in% names(priors)
  zeta_defined <- "zeta" %in% names(priors)
  
  pr_omega <- priors$omega %||% defaults$omega
  
  stopifnot_(
    !zeta_defined || pr_omega$zeta,
    c(
      "Prior for {.val zeta} cannot be used with prior on {.val omega} with
        fixed scale.",
      i = "Use, e.g., {.code omega = logistic_normal_pr(NA)}, or drop the prior
      for {.val zeta}."
    )
  )
  if (!pr_omega$zeta) defaults$zeta <- NULL
  invalid_priors <- setdiff(names(priors), names(defaults))
  stopifnot_(
    length(invalid_priors) == 0,
    c(
      "Invalid parameter {cli::qty(invalid_priors)}name{?s} in {.arg priors}:
      {.val {invalid_priors}}.",
      i = "Parameters of this model are {.val {names(defaults)}}."
    )
  )
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
    n_pr <- pr$length
    n_ref <- defaults[[nm]]$length %||% length(setup$treated)
    
    if (n_pr != n_ref) {
      stopifnot_(
        n_pr == 1,
        c(
          "Length mismatch for parameter '{nm}'.",
          i = "Got {n_pr}, expected {n_ref}.",
          i = "Only scalar priors (length = 1) can be recycled."
        )
      )
      if (n_pr == 1 && n_ref > 1) {
        pr$location <- rep(pr$location, n_ref)
        pr$scale <- rep(pr$scale, n_ref)
        pr$shape <- rep(pr$shape, n_ref)
        pr$rate <- rep(pr$rate, n_ref)
        pr$shape1 <- rep(pr$shape1, n_ref)
        pr$shape2 <- rep(pr$shape2, n_ref)
        pr$concentration <- rep(pr$concentration, n_ref)
        pr$df <- rep(pr$df, n_ref)
        pr$length <- n_ref
        priors[[nm]] <- pr
      }
    }
  }
  utils::modifyList(defaults, priors)
}
