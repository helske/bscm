#' Find kappa matching median relative ESS
#' 
#' @noRd
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
    "{.arg r_ess} should be between 1 / J = {1/J} and 1."
  )
  d <- kappa_lookup |> filter(J == .env$J)
  kappa <- stats::approx(
    d$rESS, d$kappa, xout = target, ties = "ordered",
    rule = 2
  )$y
  round(kappa, 3)
}
