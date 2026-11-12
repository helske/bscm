#' Priors of the Bayesian synthetic control model
#'
#' These functions specify the prior distributions of the parameters of the
#' Bayesian synthetic control model estimated by [bscm()]. The priors are given
#' as a list in the `priors` argument of [bscm()], with elements named by the
#' parameters listed below, e.g.,
#' `priors = list(sigma = half_normal_pr(2), weights = dirichlet_pr(1))`.
#' Parameters not in the list get their default priors. Each parameter has a
#' single prior family, whose parameters can vary across the predictors or the
#' treated units. The priors of a fitted model are returned by [get_priors()],
#' and [available_priors()] lists the prior distributions.
#'
#' @details
#' # Parameters
#'
#' Each parameter is listed below with its mathematical symbol and the name of
#' its posterior draws in [as_draws.bscmfit()] and [as.data.frame.bscmfit()].
#' The parameter names are used in the `priors` argument of [bscm()] and in the
#' output of [get_priors()], [summary.bscmfit()], and [coef.bscmfit()]. The
#' labels used by the print method are given where they differ from the
#' parameter names.
#'
#' The default weakly informative priors are scaled by the following 
#' statistics of the pre-treatment period, with the resulting values rounded 
#' to two significant digits:
#'
#' * `mean_y`: Mean of the outcome of the treated unit.
#' * `sd_e`: Standard deviation of the outcome of the treated unit, or of its
#'   difference to the mean of the donors, whichever is smaller.
#' * `sd_y`: Standard deviation of the outcome, pooled within the treated
#'   units.
#' * `sd_x`: Standard deviation of the predictor, pooled within the treated
#'   units.
#' * `J`: Number of donors.
#'
#' The statistics `mean_y` and `sd_e` are computed separately for each treated
#' unit, and `sd_x` separately for each predictor.
#' 
#' ## `intercept` (\eqn{\alpha}, draws `alpha`)
#'
#' Intercept of each treated unit. Priors with support on the real line can be
#' used. The prior is defined for the intercept of the model in which the
#' donors and all predictors are centered by their pre-treatment means (draws
#' `a`), from which \eqn{\alpha} is computed. The default
#' prior is `normal_pr(mean_y, 2 * sd_e)`.
#'
#' ## `coef` (\eqn{\beta}, draws `beta`)
#'
#' Coefficients of the predictors, common to the treated units and the donors,
#' labelled by the names of the predictors in the print output, e.g., `x` for
#' the coefficient of predictor `x`. For a predictor with a time-varying
#' coefficient, \eqn{\beta} is the average level of the coefficient in the
#' pre-treatment period, see [tv()]. Priors with support on the real line can
#' be used. The default prior is `normal_pr(0, 2 * sd_y / sd_x)`.
#'
#' ## `tv_sd` (\eqn{\kappa}, draws `kappa`)
#'
#' Standard deviations of the spline coefficients of the time-varying
#' coefficients, labelled `sd(z)` for predictor `x` in the print output, see
#' [tv()]. Priors with positive support or with support on the real line
#' (truncated to positive values) can be used. The default prior is
#' `half_normal_pr(0.5 * sd_y / sd_x)`.
#'
#' ## `tv_slope` (\eqn{\nu}, draws `nu`)
#'
#' Slopes of the linear trend of the time-varying coefficients with
#' `type = "rw2"` in [tv()], labelled `slope(z)` for predictor `x` in the
#' print output. Priors with support on the real line can be used. The default
#' prior is `normal_pr(0, 2 * sqrt(2) / (T - 1) * sd_y / sd_x)`, where `T` is 
#' the maximum length of the pre-treatment period of the treated units.
#'
#' ## `sigma` (\eqn{\sigma}, draws `sigma`)
#'
#' Standard deviation of the error term of each treated unit. Priors with
#' positive support or with support on the real line (truncated to positive
#' values) can be used. With `marginalize_sigma = TRUE` in [bscm()], only
#' `half_normal_pr()`, `scaled_inv_chi_pr()`, and `sqrt_inv_gamma_pr()` are
#' supported. The default prior is `half_normal_pr(sd_e)`.
#'
#' ## `ar1` (\eqn{\rho}, draws `rho`)
#'
#' Autoregressive coefficient of the error term of each treated unit, with
#' `error = "ar1"` in [bscm()]. The prior is defined for the transformation
#' `0.5 * (1 + rho)`, so that its support is (0, 1), and only `beta_pr()` is
#' supported. The default prior is `beta_pr(2, 2)`.
#'
#' Donor weights of each treated unit. The prior must be one of the simplex 
#' priors `dirichlet_pr()`, `logistic_normal_pr()`, or `logistic_laplace_pr()`,
#' whose arguments `concentration` and `tilt` define the concentration and the
#' tilt of the prior, see below and the section "Prior distributions". The
#' default prior is
#' `logistic_normal_pr(tilt = list(R = 3, scale = 2))`, i.e., the logistic
#' normal prior with estimated concentration (see below) and a tilt along the 
#' first three principal component directions of the donors. With
#' `marginalize_sigma = FALSE` in [bscm()], the default prior has no tilt. A
#' user-defined prior of `weights` has a tilt only if its argument `tilt` is
#' given, e.g., `logistic_normal_pr()` defines the prior without a tilt.
#'
#' ## `concentration` (\eqn{\zeta}, draws `zeta`)
#'
#' Concentration of the prior of the donor weights of each treated unit. The
#' concentration is either fixed by giving positive number(s) as the argument
#' `concentration` of the simplex prior, or estimated by giving its prior as
#' `concentration`, e.g., `logistic_normal_pr(inv_gamma_pr(10, 2))`. Priors
#' with positive support or with support on the real line (truncated to
#' positive values) can be used. Fixed and estimated concentrations cannot be
#' mixed across treated units. The default prior (`concentration = NULL`) is
#' `inv_gamma_pr(sqrt(J), 0.5 * sqrt(J) / log(1 + J / 5))`. This prior
#' penalizes very small concentrations, which are prone to cause sampling
#' problems due to challenging posterior geometry.
#'
#' ## `tilt` (\eqn{\phi}, draws `phi`)
#'
#' Tilt coefficients of the prior of the donor weights, labelled `tilt(PC1)`
#' in the print output. The argument `tilt` of the simplex priors adds an
#' exponential tilt to the prior of the donor weights. With the tilt, the
#' weights of treated unit \eqn{i} are
#' \eqn{\omega_i = \text{softmax}(V \phi_i + u_i)}, where \eqn{u_i} are the
#' log-weights of the untilted prior, \eqn{V} is a \eqn{J \times R} matrix
#' of directions shared by all treated units, and
#' \eqn{\phi_i \sim N(0, \text{scale}^2 I)} are the tilt coefficients of
#' treated unit \eqn{i}. The tilt is given as a list, either
#'
#' * `list(R = <integer>, scale = <number>)`: The columns of \eqn{V} are the
#'   first `R` principal component directions of the donors, i.e., the first
#'   `R` right singular vectors of the matrix of the donors' pre-treatment
#'   outcomes (the first `min(T_pre)` time points) after two-way demeaning,
#'   with each column of \eqn{V} standardized to zero mean and unit standard
#'   deviation. The matrix \eqn{V} is computed from the data of each fit,
#'   so refits such as those of [placebo_effects()], [leave_donor_out()],
#'   and [lfo()] compute it from their own donors and pre-treatment period.
#' * `list(V = <matrix>, scale = <number>)`: A user-defined numeric matrix
#'   \eqn{V} with one column per direction and row names containing the ids
#'   of all donors. The rows are matched to the donors by the row names, and
#'   rows of other units are ignored. The matrix is used as is, without
#'   centering or scaling.
#'
#' In both cases, `scale` is the prior standard deviation of the tilt
#' coefficients \eqn{\phi_i}. In the output, the tilt coefficients are named
#' `PC1`, `PC2`, ... for the principal component directions, and by the column
#' names of a user-defined `V` (`V1`, `V2`, ... if it has none). The tilt
#' requires `marginalize_sigma = TRUE` in [bscm()]. The default prior of
#' `weights` has the tilt `list(R = 3, scale = 2)`.
#'
#' # Prior distributions
#'
#' The following priors follow the definitions of the corresponding
#' distributions in Stan:
#'
#' * `beta_pr(shape1, shape2)`: Beta distribution.
#' * `cauchy_pr(location, scale)`: Cauchy distribution.
#' * `exponential_pr(rate)`: Exponential distribution.
#' * `gamma_pr(shape, rate)`: Gamma distribution.
#' * `inv_gamma_pr(shape, scale)`: Inverse gamma distribution.
#' * `laplace_pr(location, scale)`: Laplace (double exponential) distribution.
#' * `normal_pr(location, scale)`: Normal distribution.
#' * `student_pr(df, location, scale)`: Student's t distribution.
#'
#' In addition, the following priors are available:
#'
#' * `half_normal_pr(scale)`: Half-normal distribution. This corresponds to
#'   `normal_pr(0, scale)` prior for a parameter with lower bound of zero.
#' * `sqrt_inv_gamma_pr(shape, scale)`: Square-root inverse gamma
#'   distribution defined as
#'   \eqn{x^2 \sim \text{Inv-Gamma}(\text{shape}, \text{scale})}, with the
#'   inverse gamma distribution as in `inv_gamma_pr()`.
#' * `scaled_inv_chi_pr(df, scale)`: Scaled inverse chi distribution defined
#'   as
#'   \eqn{x^2 \sim \text{Inv-Gamma}(\text{df} / 2, \text{scale}^2 / 2)}, i.e.,
#'   `sqrt_inv_gamma_pr(df / 2, scale^2 / 2)`.
#'
#' Of these, `cauchy_pr()`, `laplace_pr()`, `normal_pr()`, and `student_pr()`
#' have support on the real line, `beta_pr()` on (0, 1), and the others
#' positive support.
#'
#' For the donor weights \eqn{\omega}, the following simplex priors are
#' available:
#'
#' * `dirichlet_pr(concentration, tilt)`: Symmetric Dirichlet distribution
#'   with concentration \eqn{\zeta}, i.e.,
#'   \eqn{\omega \sim \text{Dirichlet}(\zeta, \ldots, \zeta)}.
#' * `logistic_normal_pr(concentration, tilt)`: Logistic normal distribution
#'   where the concentration \eqn{\zeta} is the precision of the normal
#'   distribution, so that \eqn{\omega = \text{softmax}(\eta / \sqrt{\zeta})},
#'   where \eqn{\eta \sim N(0, I)} constrained to sum to zero and scaled to
#'   have marginal standard deviation of one.
#' * `logistic_laplace_pr(concentration, tilt)`: As the logistic normal prior,
#'   but with the zero-sum vector \eqn{\eta} constructed from independent
#'   (scaled) Laplace-distributed variables.
#'
#' For all three simplex priors, small concentrations correspond to a prior
#' favoring sparse weight vectors, while large concentrations correspond to
#' more diffuse weight vectors. Specifically for the Dirichlet prior,
#' \eqn{\zeta < 1} concentrates the weight on few donors, \eqn{\zeta = 1} is
#' the uniform prior over the simplex, and \eqn{\zeta > 1} pulls the weights
#' toward the center of the simplex.
#'
#' @param concentration \[`numeric()`, `bscm_prior`, or `NULL`]\cr Fixed
#'   positive concentration parameter(s), a `bscm_prior`, or `NULL` (the
#'   default) for the default prior, see details.
#' @param df \[`numeric()`]\cr Degrees of freedom parameter(s).
#' @param location \[`numeric()`]\cr Location parameter(s).
#' @param rate \[`numeric()`]\cr Rate parameter(s).
#' @param scale \[`numeric()`]\cr Scale parameter(s).
#' @param shape \[`numeric()`]\cr Shape parameter(s).
#' @param shape1 \[`numeric()`]\cr First shape parameter(s) of the beta prior.
#' @param shape2 \[`numeric()`]\cr Second shape parameter(s) of the beta
#'   prior.
#' @param tilt \[`list()` or `NULL`]\cr Tilt of the prior of the donor
#'   weights, either `list(R = <integer>, scale = <number>)` or
#'   `list(V = <matrix>, scale = <number>)`. Default is `NULL` (no tilt). See
#'   details.
#' @return A `bscm_prior` object.
#' @seealso [available_priors()], [bscm()], [get_priors()]
#' @references Stan Development Team. 2026. Stan Functions Reference, 2.40.
#' https://mc-stan.org/docs/functions-reference
#' @examples
#' normal_pr(c(0, 0), c(1, 2))
#' student_pr(df = 4, location = 0, scale = 2)
#' exponential_pr(rate = c(1, 0.1))
#' scaled_inv_chi_pr(df = 4, scale = 1)
#' beta_pr(2, 2)
#' dirichlet_pr(1)
#' logistic_normal_pr(2)
#' logistic_normal_pr()
#' logistic_normal_pr(concentration = inv_gamma_pr(10, 2))
#' logistic_normal_pr(tilt = list(R = 3, scale = 2))
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
dirichlet_pr <- function(concentration = NULL, tilt = NULL) {
  simplex_prior("dirichlet", concentration, tilt)
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
logistic_laplace_pr <- function(concentration = NULL, tilt = NULL) {
  simplex_prior("logistic_laplace", concentration, tilt)
}
#' @export
#' @rdname bscm_prior
logistic_normal_pr <- function(concentration = NULL, tilt = NULL) {
  simplex_prior("logistic_normal", concentration, tilt)
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

#' Construct a prior on the simplex for the donor weights
#'
#' @param distribution \[`character(1)`]\cr Name of the simplex prior.
#' @param concentration \[`numeric()`, `bscm_prior`, or `NULL`]\cr Fixed
#'   concentration(s), prior of the estimated concentration, or `NULL` for
#'   estimated concentration with the default prior.
#' @param tilt \[`list()` or `NULL`]\cr Specification of the tilt.
#' @param call The calling environment reported in the error message.
#' @noRd
simplex_prior <- function(
    distribution,
    concentration,
    tilt,
    call = rlang::caller_env()
) {
  if (is.null(concentration) || inherits(concentration, "bscm_prior")) {
    n <- concentration$length %||% 1L
  } else {
    stopifnot_(
      checkmate::test_numeric(
        concentration,
        lower = 0,
        any.missing = FALSE,
        min.len = 1L,
        finite = TRUE
      ) &&
        all(concentration > 0),
      c(
        "Argument {.arg concentration} must be positive number(s), a
        {.cls bscm_prior} object, or {.code NULL}.",
        i = "Use a {.cls bscm_prior} object such as
        {.code inv_gamma_pr(2, 1)} or {.code NULL} to estimate the
        concentration."
      ),
      call = call
    )
    n <- length(concentration)
  }
  structure(
    list(
      distribution = distribution,
      concentration = concentration,
      npar = 1,
      length = n,
      tilt = check_tilt(tilt, call = call)
    ),
    class = "bscm_prior"
  )
}

#' Check and normalize the tilt specification of a simplex prior
#'
#' The tilt is given either by the number of directions `R` computed from the
#' donors, or by a user-supplied matrix `V`, together with the prior standard
#' deviation `scale` of the tilt coefficients.
#'
#' @param tilt \[`list()` or `NULL`]\cr Specification of the tilt.
#' @param call The calling environment reported in the error message.
#' @return `NULL` or a list with elements `R`, `V` (`NULL` unless supplied by
#'   the user), and `scale`.
#' @noRd
check_tilt <- function(tilt, call = rlang::caller_env()) {
  if (is.null(tilt)) {
    return(NULL)
  }
  stopifnot_(
    checkmate::test_list(tilt, names = "unique"),
    "Argument {.arg tilt} must be a named list or {.code NULL}.",
    call = call
  )
  invalid_names <- setdiff(names(tilt), c("R", "V", "scale"))
  stopifnot_(
    length(invalid_names) == 0L,
    c(
      "Invalid {cli::qty(invalid_names)}element{?s} in {.arg tilt}:
      {.val {invalid_names}}.",
      i = "Valid elements are {.val R}, {.val V}, and {.val scale}."
    ),
    call = call
  )
  has_R <- !is.null(tilt$R)
  has_V <- !is.null(tilt$V)
  stopifnot_(
    xor(has_R, has_V),
    "Argument {.arg tilt} must contain either element {.val R} or element
    {.val V}, but not both.",
    call = call
  )
  stopifnot_(
    checkmate::test_number(tilt$scale, lower = 0, finite = TRUE) &&
      tilt$scale > 0,
    "Element {.val scale} of {.arg tilt} must be a single positive number.",
    call = call
  )
  if (has_R) {
    stopifnot_(
      checkmate::test_count(tilt$R, positive = TRUE),
      "Element {.val R} of {.arg tilt} must be a single positive integer.",
      call = call
    )
    R <- as.integer(tilt$R)
  } else {
    stopifnot_(
      checkmate::test_matrix(
        tilt$V,
        mode = "numeric",
        any.missing = FALSE,
        min.rows = 2L,
        min.cols = 1L
      ) &&
        all(is.finite(tilt$V)),
      "Element {.val V} of {.arg tilt} must be a numeric matrix with finite
      values.",
      call = call
    )
    stopifnot_(
      !is.null(rownames(tilt$V)) && !anyDuplicated(rownames(tilt$V)),
      c(
        "Element {.val V} of {.arg tilt} must have unique row names.",
        i = "The row names of {.val V} are matched to the donor ids."
      ),
      call = call
    )
    R <- ncol(tilt$V)
  }
  list(R = R, V = tilt$V, scale = tilt$scale)
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
    dirichlet = simplex_prior_label(x),
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
    logistic_laplace = simplex_prior_label(x),
    logistic_normal = simplex_prior_label(x),
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
#' Plain-text label of a simplex prior
#'
#' @param x \[`bscm_prior`]\cr A prior on the simplex.
#' @noRd
simplex_prior_label <- function(x) {
  # NULL is the default prior not yet resolved by define_priors()
  concentration <- x$concentration %||% "NULL"
  if (inherits(concentration, "bscm_prior")) {
    concentration <- as.character(concentration)
  }
  pars <- paste0("concentration = ", concentration)
  if (!is.null(x$tilt)) {
    pars <- paste0(
      pars,
      ", tilt = list(R = ", x$tilt$R, ", scale = ", x$tilt$scale, ")"
    )
  }
  paste0(x$distribution, "(", pars, ")")
}
#' @export
print.bscm_prior <- function(x, ...) {
  cat(paste(as.character(x), collapse = "\n"), "\n")
  invisible(x)
}

