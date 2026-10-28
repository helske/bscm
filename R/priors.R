#' Prior distributions in the Bayesian synthetic control model
#'
#' These functions specify prior distributions for parameters of the Bayesian
#' synthetic control model estimated by [bscm()]. See [available_priors()] for
#' the list of available prior distributions and their definitions.
#'
#' Priors with support on the real line can be used for parameters that can
#' take any real value (`alpha`, `beta`, `nu`). Note however that for the
#' intercept the prior corresponds to a value of the intercept when donors and
#' predictors corresponding to time-constant coefficients are centered.
#'
#' For parameters constrained to be positive (`sigma`, `kappa`, `zeta`),
#' priors with positive support can be used, as well as the priors with
#' support on the real line, which are then truncated to positive values.
#'
#' For autoregressive coefficients of the residuals (rho) the prior corresponds
#' to transformation `0.5 * (1 + rho)`, so that the support is (0, 1). Only
#' supported prior distribution for this is `beta_pr(shape1, shape2)`.
#'
#' For the donor weight vector `omega`, priors with support on the simplex are
#' used. The concentration or scale parameter \eqn{\zeta} of these priors can
#' either be fixed by the user or estimated from the data. Calling
#' `logistic_normal_pr()` (and equivalently for `dirichlet_pr()` and
#' `logistic_laplace_pr()`) without arguments (or argument value `NA`) leaves
#' the prior parameter unspecified, in which case \eqn{\zeta_i} is a parameter
#' of the model, estimated separately for each treated unit \eqn{i}, with the 
#' prior for \eqn{\zeta_i} is defined as the element
#' `zeta` of the `priors` list of [bscm()]. The default prior for \eqn{\zeta_i} 
#' is `inv_gamma_pr(sqrt(J), 0.5 * sqrt(J) / log(1 + J / 5)`, where `J` is the 
#' number of donors. This prior penalizes very small concentrations 
#' which are prone to cause sampling problems due to challenging posterior 
#' geometry and concentration of weights for a single donor.
#' Note: This default is chosen based on prior calibration 
#' of the distribution of effective number of donors \eqn{1/\sum(\omega^2)} 
#' with the total number of donors \eqn{J} ranging from 20 to 200, 
#' leading to 5% and 95% prior quantiles about 2 and J/10.
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
#' | `zeta` | `concentration` | concentration of the prior of `omega` |
#' | `sigma` | `sd` | standard deviation of the error term |
#'
#' @param concentration \[`numeric()`]\cr Concentration parameter(s). Can
#'   also be `NA`, in which case the concentration is estimated from the data,
#'   see details.
#' @param df \[`numeric()`]\cr Degrees of freedom parameter(s).
#' @param location \[`numeric()`]\cr Location parameter(s).
#' @param rate \[`numeric()`]\cr Rate parameter(s).
#' @param scale \[`numeric()`]\cr Scale parameter(s). For the priors of the
#'   donor weights `omega`, this can also be `NA`, in which case the scale is
#'   estimated from the data, see details.
#' @param shape \[`numeric()`]\cr Shape parameter(s).
#' @param shape1 \[`numeric()`]\cr First shape parameter(s) of the beta prior.
#' @param shape2 \[`numeric()`]\cr Second shape parameter(s) of the beta
#'   prior.
#'
#' @return A `bscm_prior` object.
#' @seealso [available_priors()], [bscm()]
#' @examples
#' normal_pr(c(0, 0), c(1, 2))
#' student_pr(df = 4, location = 0, scale = 2)
#' exponential_pr(rate = c(1, 0.1))
#' scaled_inv_chi_pr(df = 4, scale = 1)
#' beta_pr(2, 2)
#' dirichlet_pr(1)
#' logistic_normal_pr(2)
#' logistic_normal_pr()
#'
#' @name bscm_prior
NULL

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
dirichlet_pr <- function(concentration = NA) {
  zeta <- anyNA(concentration)
  if (zeta) {
    stopifnot_(
      all(is.na(concentration)),
      "Argument {.arg concentration} of {.fun dirichlet_pr} must be 
      either all positive values or all {.code NA}, in which case the 
      concentrations are defined their own priors using the element named 
      {.val zeta} of the list argument {.arg prior}."
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
inv_gamma_pr <- function(shape, scale) {
  check_positive(shape, "shape")
  check_positive(scale, "scale")
  n <- max(length(shape), length(scale))
  if (length(shape) == 1) {
    shape <- rep_len(shape, n)
  }
  if (length(scale) == 1) {
    scale <- rep_len(scale, n)
  }
  stopifnot_(
    length(shape) == length(scale),
    "Lengths of {.arg shape} and {.arg scale} should match."
  )
  structure(
    list(
      distribution = "inv_gamma",
      shape = shape,
      scale = scale,
      npar = 2,
      length = n
    ),
    class = "bscm_prior"
  )
}
#' @export
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
logistic_laplace_pr <- function(concentration = NA) {
  zeta <- anyNA(concentration)
  if (zeta) {
    stopifnot_(
      all(is.na(concentration)),
      "Argument {.arg concentration} of {.fun logistic_laplace_pr} must be 
      either all positive values or all {.code NA}, in which case the 
      concentrations are defined their own priors using the element named 
      {.val zeta} of the list argument {.arg prior}."
    )
    concentration <- rep(NA, length(concentration))
  } else {
    check_positive(concentration, "concentration")
  }
  structure(
    list(
      distribution = "logistic_laplace",
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
logistic_normal_pr <- function(concentration = NA) {
  zeta <- anyNA(concentration)
  if (zeta) {
    stopifnot_(
      all(is.na(concentration)),
      "Argument {.arg concentration} of {.fun logistic_normal_pr} must be 
      either all positive values or all {.code NA}, in which case the 
      concentrations are defined their own priors using the element named 
      {.val zeta} of the list argument {.arg prior}."
    )
    concentration <- rep(NA, length(concentration))
  } else {
    check_positive(concentration, "concentration")
  }
  structure(
    list(
      distribution = "logistic_normal",
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
#' @export
#' @rdname bscm_prior
scaled_inv_chi_pr <- function(df, scale) {
  check_positive(df, "df")
  check_positive(scale, "scale")
  n <- max(length(df), length(scale))
  if (length(df) == 1) {
    df <- rep_len(df, n)
  }
  if (length(scale) == 1) {
    scale <- rep_len(scale, n)
  }
  stopifnot_(
    length(df) == length(scale),
    "Lengths of {.arg df} and {.arg scale} should match."
  )
  structure(
    list(
      distribution = "scaled_inv_chi",
      df = df,
      scale = scale,
      npar = 2,
      length = n
    ),
    class = "bscm_prior"
  )
}
#' @export
#' @rdname bscm_prior
sqrt_inv_gamma_pr <- function(shape, scale) {
  check_positive(shape, "shape")
  check_positive(scale, "scale")
  n <- max(length(shape), length(scale))
  if (length(shape) == 1) {
    shape <- rep_len(shape, n)
  }
  if (length(scale) == 1) {
    scale <- rep_len(scale, n)
  }
  stopifnot_(
    length(shape) == length(scale),
    "Lengths of {.arg shape} and {.arg scale} should match."
  )
  structure(
    list(
      distribution = "sqrt_inv_gamma",
      shape = shape,
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
  check_positive(df, "df")
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
as.character.bscm_prior <- function(x, ...) {

  switch(
    x$distribution,
    beta = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0("beta(", x$shape1[i], ", ", x$shape2[i], ")")
        },
        character(1)
      )
    },
    cauchy = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0("cauchy(", x$location[i], ", ", x$scale[i], ")")
        },
        character(1)
      )
    },
    dirichlet = {
      if (x$zeta) {
        "dirichlet(zeta)"
      } else {
        paste0("dirichlet(", x$concentration, ")")
      }
    },
    exponential = {
      paste0("exponential(", x$rate, ")")
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
    half_normal = {
      paste0("half_normal(", x$scale, ")")
    },
    inv_gamma = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0("inv_gamma(", x$shape[i], ", ", x$scale[i], ")")
        },
        character(1)
      )
    },
    laplace = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0("laplace(", x$location[i], ", ", x$scale[i], ")")
        },
        character(1)
      )
    },
    logistic_laplace = {
      if (x$zeta) {
        "logistic_laplace(zeta)"
      } else {
        paste0("logistic_laplace(", x$concentration, ")")
      }
    },
    logistic_normal = {
      if (x$zeta) {
        "logistic_normal(zeta)"
      } else {
        paste0("logistic_normal(", x$concentration, ")")
      }
    },
    normal = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0("normal(", x$location[i], ", ", x$scale[i], ")")
        },
        character(1)
      )
    },
    scaled_inv_chi = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0("scaled_inv_chi(", x$df[i], ", ", x$scale[i], ")")
        },
        character(1)
      )
    },
    sqrt_inv_gamma = {
      vapply(
        seq_len(x$length),
        function(i) {
          paste0("sqrt_inv_gamma(", x$shape[i], ", ", x$scale[i], ")")
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
    }
  )
}

#' @export
print.bscm_prior <- function(x, ...) {
  cat(paste(as.character(x), collapse = "\n"), "\n")
  invisible(x)
}
