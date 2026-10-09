#' Prior distributions in the Bayesian synthetic control model
#'
#' The functions `normal_pr()`, `student_pr()`, `exponential_pr()`, `gamma_pr()`,
#' `beta_pr()`, `half_normal_pr()`, `dirichlet_pr()`, and `logistic_normal_pr()`
#' specify prior distributions for parameters of the Bayesian synthetic
#' control model estimated by [bscm()].
#'
#' For parameters that can take any real value (`alpha`, `beta`, `nu`),
#' currently supported prior distributions are `normal_pr(location, scale)` and
#' `student_pr(df, location, scale)`. Note however that for the intercept the
#' prior corresponds to a value of the intercept when donors and predictors
#' corresponding to time-constant coefficients are centered.
#'
#' For parameters constrained to be positive (`sigma`, `kappa`, `zeta`),
#' supported distributions are `exponential_pr(rate)`, `gamma_pr(shape, rate)`,
#' and `half_normal_pr(scale)`.
#'
#' For autoregressive coefficients of the residuals (rho) the prior corresponds
#' to transformation `0.5 * (1 + rho)`, so that the support is (0, 1). Only
#' supported prior distribution for this is `beta_pr(shape1, shape2)`.
#'
#' For the donor weight vector `omega`, supported priors are
#' `dirichlet_pr(concentration)` and `logistic_normal_pr(scale)`. The
#' (symmetric) Dirichlet prior with concentration \eqn{\zeta} defined as
#' \eqn{\omega \sim \text{Dirichlet}(\zeta, \ldots,\zeta)}.
#' Values \eqn{\zeta < 1} concentrate weight on few donors, while
#' \eqn{\zeta = 1} is uniform over the simplex, and
#' \eqn{\zeta > 1} pulls weights toward the center of the simplex.
#' Logistic normal with with scale \eqn{\zeta} is defined as
#' \eqn{\omega = \text{softmax}(\zeta\eta)},
#' where \eqn{\eta \sim N(0, I)} constrained to sum to zero. Larger \eqn{\zeta}
#' induces more concentrated (sparser) weights, while \eqn{\zeta \to 0}
#' corresponds to the uniform weights \eqn{1/J}.
#'
#' The \eqn{\zeta} parameter of the Dirichlet and logistic normal priors can 
#' either be fixed by the user or estimated from the data. Calling 
#' `logistic_normal_pr()` (and equivalently for `dirichlet_pr()`without 
#' arguments (or argument value `NA`) leaves the prior parameter unspecified, in
#' which case \eqn{\zeta_i} is a parameter of the model, estimated separately
#' for each treated unit \eqn{i}. When the scale is estimated, the prior for 
#' \eqn{\zeta_i} is defined as the element `zeta` of the `priors` list of 
#' [bscm()].
#'
#' Each element of the `priors` list of [bscm()] can be named either by the
#' name of the model parameter or by a plain-text synonym:
#'
#' | Parameter | Synonym | Quantity |
#' | --- | --- | --- |
#' | `alpha` | `intercept` | intercept of each treated unit |
#' | `beta` | `coef` | time-constant regression coefficients |
#' | `kappa` | `scale_gamma` | scale of the time-varying coefficients |
#' | `nu` | `slope_gamma` | trend slope of the time-varying coefficients |
#' | `rho` | `ar1` | autoregressive coefficient of the residuals |
#' | `omega` | `weights` | donor weights |
#' | `zeta` | `scale_omega` | scale or concentration of Dirirchlet or logistic normal prior |
#' | `sigma` | `sd` | standard deviation of the error term |
#'
#'
#' @param location \[`numeric(1)`]\cr Location parameter for normal and
#'   Student-t priors.
#' @param scale \[`numeric(1)`]\cr Scale parameter for normal, Student-t,
#'   half-normal, and logistic normal priors. For the logistic normal prior
#'   this can also be `NA`, in which case the scale is estimated
#'   from the data, see details.
#' @param df \[`numeric(1)`]\cr Degrees of freedom for the Student-t prior.
#' @param shape \[`numeric(1)`]\cr Shape parameter for the gamma prior.
#' @param rate \[`numeric(1)`]\cr Rate parameter for the gamma and exponential
#'   priors.
#' @param shape1 \[`numeric(1)`]\cr First shape parameter for the beta prior.
#' @param shape2 \[`numeric(1)`]\cr Second shape parameter for the beta prior.
#' @param concentration \[`numeric(1)`]\cr Concentration parameter for the
#'   Dirichlet prior.
#'
#' @return A `bscm_prior` object.
#' @seealso [bscm()]
#' @examples
#' normal_pr(c(0, 0), c(1, 2))
#' student_pr(df = 4, location = 0, scale = 2)
#' exponential_pr(rate = c(1, 0.1))
#' gamma_pr(2, 1)
#' beta_pr(2, 2)
#' half_normal_pr(1)
#' dirichlet_pr(1)
#' logistic_normal_pr(2)
#' logistic_normal_pr()
#'
#' @name bscm_prior
NULL
#' @param location The location parameter for the prior.
#' @param scale The scale parameter for the prior.
#' @param df The degrees of freedom for the Student's t prior.
#' @param shape The shape parameter for the gamma prior.
#' @param rate The rate parameter for the gamma and exponential priors.
#' @param shape1 The first shape parameter for the beta prior.
#' @param shape2 The second shape parameter for the beta prior.
#' @param concentration The concentration parameter for the Dirichlet prior.
#' @export
#' @rdname bscm_prior
beta_pr <- function(shape1, shape2) {
  check_positive(shape1, "shape1")
  check_positive(shape2, "shape2")
  n <- max(length(shape1), length(shape2))
  if (length(shape1) == 1) {
    shape1 <- rep_len(shape1, n)
  }
  if (length(shape2) == 1) {
    shape2 <- rep_len(shape2, n)
  }
  stopifnot_(
    length(shape1) == length(shape2),
    "Lengths of {.arg shape1} and {.arg shape2} should match."
  )
  structure(
    list(
      distribution = "beta",
      shape1 = shape1,
      shape2 = shape2,
      npar = 2,
      length = n
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
cauchy_pr <- function(location, scale) {
  check_real(location, "location")
  check_positive(scale, "scale")
  n <- max(length(location), length(scale))
  if (length(location) == 1) {
    location <- rep_len(location, n)
  }
  if (length(scale) == 1) {
    scale <- rep_len(scale, n)
  }
  stopifnot_(
    length(location) == length(scale),
    "Lengths of {.arg location} and {.arg scale} should match."
  )
  structure(
    list(
      distribution = "cauchy",
      location = location,
      scale = scale,
      npar = 2,
      length = n
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
exponential_pr <- function(rate) {
  check_positive(rate, "rate")
  structure(
    list(
      distribution = "exponential",
      rate = rate,
      npar = 1,
      length = length(rate)
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
gamma_pr <- function(shape, rate) {
  check_positive(shape, "shape")
  check_positive(rate, "rate")
  n <- max(length(shape), length(rate))
  if (length(shape) == 1) {
    shape <- rep_len(shape, n)
  }
  if (length(rate) == 1) {
    rate <- rep_len(rate, n)
  }
  stopifnot_(
    length(shape) == length(rate),
    "Lengths of {.arg shape} and {.arg rate} should match."
  )
  structure(
    list(
      distribution = "gamma",
      shape = shape,
      rate = rate,
      npar = 2,
      length = n
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
half_normal_pr <- function(scale) {
  check_positive(scale, "scale")
  
  structure(
    list(
      distribution = "half_normal",
      scale = scale,
      npar = 1,
      length = length(scale)
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
normal_pr <- function(location, scale) {
  check_real(location, "location")
  check_positive(scale, "scale")
  n <- max(length(location), length(scale))
  if (length(location) == 1) {
    location <- rep_len(location, n)
  }
  if (length(scale) == 1) {
    scale <- rep_len(scale, n)
  }
  stopifnot_(
    length(location) == length(scale),
    "Lengths of {.arg location} and {.arg scale} should match."
  )
  structure(
    list(
      distribution = "normal",
      location = location,
      scale = scale,
      npar = 2,
      length = n
    ),
    class = "bscm_prior"
  )
}
#' @rdname bscm_prior
laplace_pr <- function(location, scale) {
  check_real(location, "location")
  check_positive(scale, "scale")
  n <- max(length(location), length(scale))
  if (length(location) == 1) {
    location <- rep_len(location, n)
  }
  if (length(scale) == 1) {
    scale <- rep_len(scale, n)
  }
  stopifnot_(
    length(location) == length(scale),
    "Lengths of {.arg location} and {.arg scale} should match."
  )
  structure(
    list(
      distribution = "laplace",
      location = location,
      scale = scale,
      npar = 2,
      length = n
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
student_pr <- function(df, location, scale) {
  check_positive(scale, "df")
  check_real(location, "location")
  check_positive(scale, "scale")
  n <- max(length(df), length(location), length(scale))
  if (length(df) == 1) {
    df <- rep_len(df, n)
  }
  if (length(location) == 1) {
    location <- rep_len(location, n)
  }
  if (length(scale) == 1) {
    scale <- rep_len(scale, n)
  }
  stopifnot_(
    length(location) == length(scale) && length(location) == length(df),
    "Lengths of {.arg location}, {.arg scale}, and {.arg df} should match."
  )
  structure(
    list(
      distribution = "student_t",
      df = df,
      location = location,
      scale = scale,
      npar = 3,
      length = n
    ),
    class = "bscm_prior"
  )
}

#' @export
#' @rdname bscm_prior
dirichlet_pr <- function(concentration = NA) {
  zeta <- anyNA(concentration)
  if (zeta) {
    stopifnot_(
      all(is.na(concentration)),
      "Argument {.arg concentration} of {.fun dirichlet_pr} must be either all 
      positive values or all {.code NA}, in which case the concentrations are 
      defined priors using the element named {.val zeta} of the list argument 
      {.arg prior}."
    )
    concentration <- rep(NA, length(concentration))
  } else {
    check_positive(concentration, "concentration")
  }
  
  structure(
    list(
      distribution = "dirichlet",
      concentration = concentration,
      npar = 1,
      length = length(concentration),
      zeta = zeta
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
logistic_t_pr <- function(scale = NA) {
  # Note: The degrees of freedom is hard-coded to 3
  zeta <- anyNA(scale)
  if (zeta) {
    stopifnot_(
      all(is.na(scale)),
      "Argument {.arg scale} of {.fun logistic_t_pr} must be either all 
      positive values or all {.code NA}, in which case the scales are 
      defined priors using the element named {.val zeta} of the list argument 
      {.arg prior}."
    )
    scale <- rep(NA, length(scale))
  } else {
    check_positive(scale, "scale")
  }
  structure(
    list(
      distribution = "logistic_t",
      scale = scale,
      npar = 1,
      length = length(scale),
      zeta = zeta
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
logistic_laplace_pr <- function(scale = NA) {
  zeta <- anyNA(scale)
  if (zeta) {
    stopifnot_(
      all(is.na(scale)),
      "Argument {.arg scale} of {.fun logistic_laplace_pr} must be either all 
      positive values or all {.code NA}, in which case the scales are 
      defined priors using the element named {.val zeta} of the list argument 
      {.arg prior}."
    )
    scale <- rep(NA, length(scale))
  } else {
    check_positive(scale, "scale")
  }
  structure(
    list(
      distribution = "logistic_laplace",
      scale = scale,
      npar = 1,
      length = length(scale),
      zeta = zeta
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
logistic_normal_pr <- function(scale = NA) {
  zeta <- anyNA(scale)
  if (zeta) {
    stopifnot_(
      all(is.na(scale)),
      "Argument {.arg scale} of {.fun logistic_normal_pr} must be either all 
      positive values or all {.code NA}, in which case the scales are 
      defined priors using the element named {.val zeta} of the list argument 
      {.arg prior}."
    )
    scale <- rep(NA, length(scale))
  } else {
    check_positive(scale, "scale")
  }
  structure(
    list(
      distribution = "logistic_normal",
      scale = scale,
      npar = 1,
      length = length(scale),
      zeta = zeta
    ),
    class = "bscm_prior"
  )
}

#' @export
as.character.bscm_prior <- function(x, ...) {
  
  switch(
    x$distribution,
    normal = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0("normal(", x$location[i], ", ", x$scale[i], ")")
        },
        character(1)
      )
    },
    student_t = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0(
            "student(",
            x$df[i],
            ", ",
            x$location[i],
            ", ",
            x$scale[i],
            ")"
          )
        },
        character(1)
      )
    },
    gamma = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0("gamma(", x$shape[i], ", ", x$rate[i], ")")
        },
        character(1)
      )
    },
    beta = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0("beta(", x$shape1[i], ", ", x$shape2[i], ")")
        },
        character(1)
      )
    },
    exponential = {
      paste0("exponential(", x$rate, ")")
    },
    half_normal = {
      paste0("half_normal(", x$scale, ")")
    },
    dirichlet = {
      if (x$zeta) {
        "dirichlet(zeta)"
      } else {
        paste0("dirichlet(", x$concentration, ")")
      }
    },
    logistic_normal = {
      if (x$zeta) {
        "logistic_normal(zeta)"
      } else {
        paste0("logistic_normal(", x$scale, ")")
      }
    },
    logistic_t = {
      if (x$zeta) {
        paste0("logistic_t(", x$df, ", zeta)")
      } else {
        paste0("logistic_t(", x$df, ", ", x$scale, ")")
      }
    },
    logistic_laplace = {
      if (x$zeta) {
        "logistic_laplace(zeta)"
      } else {
        paste0("logistic_laplace(", x$scale, ")")
      }
    },
  )
}

#' @export
print.bscm_prior <- function(x, ...) {
  cat(paste(as.character(x), collapse = "\n"), "\n")
  invisible(x)
}

#' @export
#' @rdname get_priors
get_priors <- function(x, ...) {
  UseMethod("get_priors")
}
#' Get prior specifications used in the Bayesian synthetic control model
#' @param x \[`bscmfit`]\cr Output from [bscm()].
#' @param ... Ignored.
#' @return A named list of `bscm_prior` objects for the parameters of the model.
#' @export
#' @rdname get_priors
get_priors.bscmfit <- function(x, ...) {
  x$priors
}

default_priors <- function(d, setup, spline_def) {
  N <- length(setup$treated)
  pr_sigma <- exponential_pr(signif(1 / d$sd_e, 2))
  pr_omega <- dirichlet_pr(rep(1, N))
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
  scale_omega = "zeta",
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
      {?is/are} given more than once in {.arg priors}.",
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
    "dirichlet", "logistic_normal", "logistic_laplace", "logistic_t"
  )
  bounded_priors <- c(
    "exponential", "gamma", "half_normal"
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
  stopifnot_(
    checkmate::test_list(priors, types = "bscm_prior", names = "named"),
    "Argument {.arg priors} must be a named list of {.cls bscm_prior} objects."
  )
  priors <- resolve_prior_aliases(priors)
  has_zeta <- "zeta" %in% names(priors)
  has_omega <- "omega" %in% names(priors)
  
  stopifnot_(
    !has_zeta || has_omega && isTRUE(priors$omega$zeta),
    c(
      "Prior for {.val zeta} cannot be used with prior on {.val omega} with 
        fixed scale.",
      i = "Use, e.g., {.code omega = logistic_normal_pr(NA)}, or drop the prior 
      for {.val zeta}."
    )
  )
  stopifnot_(
    has_zeta || !isTRUE(priors$omega$zeta),
    c(
      "Prior of {.val omega} is specified without a fixed 
      hyperparameter or a prior for {.val zeta}.",
      i = "Either fix the concentration/scale or 
      use {.code zeta} to define the prior for {.val zeta}."
    )
  )
  
  used_parameters <- c(names(defaults), if (has_zeta) zeta = "zeta")
  invalid_priors <- setdiff(names(priors), used_parameters)
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
