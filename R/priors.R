#' Prior distributions in the Bayesian synthetic control model
#'
#' The functions `normal_pr()`, `student_pr()`, `exponential_pr()`, `gamma_pr()`,
#' `beta_pr()`, `half_normal_pr()`, `dirichlet_pr()`, and `logistic_normal_pr()`
#' specify prior distributions for parameters of the Bayesian synthetic
#' control model estimated by [bscm()].
#'
#' For parameters that can take any real value (`alpha`, `beta`),
#' currently supported prior distributions are `normal_pr(location, scale)` and
#' `student_pr(df, location, scale)`. Note however that for the intercept the
#' prior corresponds to a value of the intercept when donors and predictors
#' corresponding to time-constant coefficients are centered.
#'
#' For parameters constrained to be positive (`sigma`, `sigma_gamma`),
#' supported distributions are `exponential_pr(rate)`, `gamma_pr(shape, rate)`,
#' and `half_normal_pr(scale)`.
#'
#' For autoregressive coefficients of the residuals (rho) the prior corresponds
#' to transformation `0.5 * (1 + rho)`, so that the support is (0, 1). Only
#' supported prior distribution for this is `beta_pr(shape1, shape2)`.
#'
#' For the donor weight vector `omega`, supported priors are
#' `dirichlet_pr(concentration)` and `logistic_normal_pr(scale)`. The
#' (symmetric) Dirichlet prior with concentration \eqn{\alpha} is defined as
#' \eqn{\omega \sim \text{Dirichlet}(\kappa, \ldots, \kappa)}.
#' Values \eqn{\alpha < 1} concentrate weight on few donors;
#' \eqn{\alpha = 1} is uniform over the simplex; \eqn{\alpha > 1}
#' pulls weights toward the center of the simplex.
#' Logistic normal with with scale \eqn{\sigma} is defined as
#' \eqn{\omega = \text{softmax}(\eta)},
#' where \eqn{\eta \sim N(0, \sigma^2 I)} constrained to sum to zero. Larger
#' `scale` induces more concentrated (sparser) weights.
#'
#' @param location \[`numeric(1)`]\cr Location parameter for normal and
#'   Student-t priors.
#' @param scale \[`numeric(1)`]\cr Scale parameter for normal, Student-t,
#'   half-normal, and logistic normal priors.
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
    "Lengths of {.arg location} and {.arg scale} should be match."
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
    "Lengths of {.arg location}, {.arg scale}, and {.arg df} should be match."
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
    "Lengths of {.arg shape} and {.arg rate} should be match."
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
    "Lengths of {.arg shape1} and {.arg shape2} should be match."
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
dirichlet_pr <- function(concentration) {
  check_positive(concentration, "concentration")
  
  structure(
    list(
      distribution = "dirichlet",
      concentration = concentration,
      npar = 1,
      length = length(concentration)
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
logistic_normal_pr <- function(scale) {
  check_positive(scale, "scale")
  
  structure(
    list(
      distribution = "logistic_normal",
      scale = scale,
      npar = 1,
      length = length(scale)
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
      paste0("dirichlet(", x$concentration, ")")
    },
    logistic_normal = {
      paste0("logistic_normal(", x$scale, ")")
    }
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

default_priors <- function(descriptives, setup, spline_def) {
  
  N <- length(setup$treated)
  pr_sigma <- exponential_pr(signif(1 / descriptives$sd_e, 2))
  pr_omega <- dirichlet_pr(rep(1, N))
  pr_alpha <- pr_beta <- pr_sigma_gamma <- pr_rho <- NULL
  if (setup$has_icpt) {
    pr_alpha <- normal_pr(
      signif(descriptives$mean_y, 2), 
      signif(2 * descriptives$sd_e, 2)
    )
  }
  if (setup$has_x) {
    sd_ex <- descriptives$md_sd_e / descriptives$md_sd_x
    pr_beta <- normal_pr(0, signif(2 * sd_ex, 2))
    if (setup$has_w) {
      pr_sigma_gamma <- half_normal_pr(
        signif(0.5 * sd_ex[setup$tv_idx] * spline_def$scale, 2)
      )
    }
  }
  if (setup$has_ar1) {
    pr_rho <- beta_pr(rep(2, N), rep(2, N))
  }
  priors <- list(
    sigma = pr_sigma,
    omega = pr_omega,
    intercept = pr_alpha,
    beta = pr_beta,
    sigma_gamma = pr_sigma_gamma,
    rho = pr_rho
  )
  priors[lengths(priors) > 0]
}

define_priors <- function(priors, descriptives, setup, spline_def) {
  defaults <- default_priors(descriptives, setup, spline_def)
  if (is.null(priors)) {
    return(defaults)
  }
  valid_priors <- list(
    sigma = c("exponential", "gamma", "half_normal"),
    omega = c("dirichlet", "logistic_normal"),
    intercept = c("normal", "student_t"),
    beta = c("normal", "student_t"),
    sigma_gamma = c("exponential", "gamma", "half_normal"),
    rho = "beta"
  )
  stopifnot_(
    checkmate::test_list(priors, types = "bscm_prior", names = "named"),
    "Argument {.arg priors} must be a named list of {.cls bscm_prior} objects."
  )
  invalid_priors <- setdiff(names(priors), names(defaults))
  stopifnot_(
    length(invalid_priors) == 0,
    c(
      "Invalid parameter names in {.arg priors}.",
      i = "Model parameter are {names(valid_priors)}."
    )
  )
  for (nm in names(priors)) {
    pr <- priors[[nm]]
    ref <- defaults[[nm]]
    dist <- pr$distribution
    stopifnot_(
      dist %in% valid_priors[[nm]],
      c(
        "Invalid prior family '{dist}' for parameter '{nm}'.",
        i = "Allowed families: {valid_priors[[nm]]}"
      )
    )
    n_pr <- pr$length
    n_ref <- ref$length
    
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
        pr$length <- n_ref
        priors[[nm]] <- pr
      }
    }
  }
  utils::modifyList(defaults, priors)
}
