#' Prior distribution for the donor weight vector
#'
#' `logistic_normal()` and `dirichlet()` specify the prior on the
#' simplex-valued donor weight vector \eqn{\omega} for use in [bscm()].
#'
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
#' @param kappa \[`numeric(1)`]\cr Positive concentration parameter of
#' symmetric Dirichlet distribution or positive scale parameter of logistic
#' normal distribution.
#' @return An `omega_prior` object (a list with elements `distribution` and
#'   `kappa`).
#' @seealso [bscm()]
#' @rdname omega_prior
#' @examples
#' logistic_normal(kappa = 2)
#' dirichlet(kappa = 0.5)
#' get_omega_prior(fit_single_treated)
#' @export
logistic_normal <- function(kappa) {
  check_omega_prior_args(kappa)
  structure(
    list(distribution = "logistic_normal", kappa = kappa),
    class = "omega_prior"
  )
}

#' @rdname omega_prior
#' @export
dirichlet <- function(kappa) {
  check_omega_prior_args(kappa)
  structure(
    list(distribution = "dirichlet", kappa = kappa),
    class = "omega_prior"
  )
}

#' Function `get_omega_prior()` returns the prior of weight vector \eqn{\omega}
#' used in the model estimation.
#'
#' @param x \[`bscmfit`]\cr The output returned by the [bscm()].
#' @rdname omega_prior
#' @export
get_omega_prior <- function(x) {
  stopifnot_(
    inherits(x, "bscmfit"),
    "Argument {.arg x} should be an object of class {.cls bscmfit}."
  )
  x$setup$omega_prior
}
#' @export
print.omega_prior <- function(x, ...) {
  cat(paste0(x$distribution, "(kappa = ", x$kappa, ")"), "\n")
  invisible(x)
}
#' @export
as.character.omega_prior <- function(x, ...) {
  paste0(x$distribution, "(", x$kappa, ")")
}

check_omega_prior_args <- function(kappa) {
  stopifnot_(
    checkmate::test_number(kappa, lower = 0, finite = TRUE),
    "Argument {.arg kappa} must be a single positive number."
  )
}
