#' Prior distribution for the donor weight vector
#'
#' `logistic_normal()` and `dirichlet()` specify the prior on the
#' simplex-valued donor weight vector \eqn{\omega} for use in [bscm()].
#'
#' @details
#' For the **logistic normal** prior, \eqn{\omega = \text{softmax}(\eta)}
#' where \eqn{\eta \sim N(0, \kappa^2 I)} constrained to sum to zero. Larger
#' \eqn{\kappa} induces more concentrated (sparser) weights.
#'
#' For the **symmetric Dirichlet** prior,
#' \eqn{\omega \sim \text{Dirichlet}(\kappa, \ldots, \kappa)}.
#' Values \eqn{\kappa < 1} concentrate weight on few donors;
#' \eqn{\kappa = 1} is uniform over the simplex; \eqn{\kappa > 1}
#' pulls weights toward the center of the simplex.
#'
#' Instead of supplying \eqn{\kappa} directly, a target prior median relative
#' effective sample size \eqn{r_{ESS} \in (0, 1)} can be supplied via
#' `r_ess`. Kappa selection via `r_ess` for the Dirichlet prior is not yet
#' implemented; provide `kappa` directly.
#'
#' @param kappa \[`numeric(1)`]\cr Positive concentration parameter
#'   \eqn{\kappa}. Cannot be combined with `r_ess`.
#' @param r_ess \[`numeric(1)`]\cr Target prior median relative effective
#'   number of donors, a number strictly between 0 and 1. Cannot be combined
#'   with `kappa`.
#' @return An `omega_prior` object (a list with fields `distribution`,
#'   `kappa`, and `r_ess`).
#' @seealso [bscm()]
#' @examples
#' logistic_normal(kappa = 2)
#' logistic_normal(r_ess = 0.1)
#' dirichlet(kappa = 0.5)
#' @export
logistic_normal <- function(kappa = NULL, r_ess = NULL) {
  check_omega_prior_args(kappa, r_ess)
  structure(
    list(distribution = "logistic_normal", kappa = kappa, r_ess = r_ess),
    class = "omega_prior"
  )
}

#' @rdname logistic_normal
#' @export
dirichlet <- function(kappa = NULL, r_ess = NULL) {
  check_omega_prior_args(kappa, r_ess)
  structure(
    list(distribution = "dirichlet", kappa = kappa, r_ess = r_ess),
    class = "omega_prior"
  )
}

#' @export
format.omega_prior <- function(x, ...) {
  if (!is.null(x$kappa)) {
    paste0(x$distribution, "(kappa = ", x$kappa, ")")
  } else if (!is.null(x$r_ess)) {
    paste0(x$distribution, "(r_ess = ", x$r_ess, ")")
  } else {
    paste0(x$distribution, "()")
  }
}

#' @export
print.omega_prior <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

check_omega_prior_args <- function(kappa, r_ess) {
  stopifnot_(
    checkmate::test_number(kappa, lower = 0, finite = TRUE, null.ok = TRUE),
    "Argument {.arg kappa} must be a single positive number."
  )
  stopifnot_(
    checkmate::test_number(
      r_ess,
      lower = 0,
      upper = 1,
      finite = TRUE,
      null.ok = TRUE
    ),
    "Argument {.arg r_ess} must be a single number strictly between 0 and 1."
  )
  stopifnot_(
    is.null(kappa) || is.null(r_ess),
    "Only one of {.arg kappa} and {.arg r_ess} can be specified."
  )
}

# Resolves the numeric kappa for the omega prior given J (number of donors).
resolve_kappa <- function(prior, J) {
  if (!is.null(prior$kappa)) {
    return(prior$kappa)
  }
  select_kappa(J, prior$r_ess %||% 0.25, prior$distribution)
}
