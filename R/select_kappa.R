#' Find a value for kappa by matching it to target median relative ESS
#'
#' Given a number of donors \eqn{J} and a target value for relative effective number
#' of donors, \eqn{rESS = 1 / \sum(\omega^2)}, function `select_kappa()` returns
#' a value of \eqn{\kappa} such that the prior median of  rESS given
#' \eqn{\kappa} matches the target rESS. This is based on on a lookup table
#' from Monte Carlo simulations.
#'
#' For reverse, function `get_ess()` returns the prior ESS (not relative)
#' value based on \eqn{J} and \eqn{\kappa}.
#'
#' Function `get_kappa()` returns the value of \eqn{\kappa} used in the model
#' estimation.
#'
#' @param J \[`integer(1)`]\cr Number of donors.
#' @param target \[`numeric(1)`]\cr Target value of relative ESS.
#' @param kappa \[`numeric(1)`]\cr Value of \eqn{\kappa}.
#' @param x \[`bscmfit`]\cr The output returned by the [bscm()].
#' @return For `select_kappa` and `get_kappa`, a value of \eqn{\kappa}.
#' For `get_ess()`, value of relative ESS.
#' @rdname kappa
#' @export
select_kappa <- function(J, target) {
  min_J <- min(kappa_lookup$J)
  max_J <- max(kappa_lookup$J)
  stopifnot_(
    J >= min_J && J <= max_J,
    "Defining {.arg kappa} using relative ESS via argument {.arg r_ess} is only 
    supported when the number of donors is between {min_J} and {max_J}."
  )
  stopifnot_(
    target > 1 / J && target < 1,
    c(
      "{.arg r_ess} should be between 1 / J = {round(1 / J, 4)} and 1.",
      i = "Input value of {.arg r_ess} was {target}."
    )
  )
  d <- kappa_lookup |> filter(J == .env$J)
  kappa <- stats::approx(
    d$rESS,
    d$kappa,
    xout = target,
    ties = "ordered",
    rule = 2
  )$y
  round(kappa, 3)
}
#' @rdname kappa
#' @export
get_kappa <- function(x) {
  stopifnot_(
    inherits(x, "bscmfit"),
    "Argument {.arg x} should be an object of class {.cls bscmfit}."
  )
  x$setup$kappa
}
#' @rdname kappa
#' @export
get_ess <- function(J, kappa) {
  min_J <- min(kappa_lookup$J)
  max_J <- max(kappa_lookup$J)
  stopifnot_(
    J >= min_J && J <= max_J,
    "Extracting prior relative ESS by {.arg kappa} is only 
    supported when the number of donors is between {min_J} and {max_J}."
  )
  stopifnot_(
    kappa >= 0.05 && kappa <= 10,
    "Extracting prior ESS by {.arg kappa} is only 
    supported when {.arg kappa} is between 0.05 and 10."
  )
  d <- kappa_lookup |> filter(J == .env$J) |> arrange(kappa)
  r_ess <- stats::approx(
    d$kappa,
    d$rESS,
    xout = kappa,
    ties = "ordered",
    rule = 2
  )$y
  round(J * r_ess, 2)
}
