#' @export
#' @rdname leave_donor_out
leave_donor_out <- function(x, ...) {
  UseMethod("leave_donor_out", x)
}
#' Leave-donor-out effects of a Bayesian synthetic control model
#' 
#' Re-estimates original model \eqn{J} times by removing one of the \eqn{J} 
#' donors from the data in turn. This can be used to assess how sensitive the 
#' results are to the inclusion of specific donors.
#'   
#' @export
#' @rdname leave_donor_out
#' @param x \[`bscmfit`]\cr The output returned by the [bscm()].
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries of the 
#' treatment effects and RMSE estimates. Default is `c(0.025, 0.975)`.
#' @param ... Additional arguments passed on to [bscm()].
#' @return A list with three elements: `effect`, `rmse`, and `diagnostics`, 
#' containing the treatment effect estimates, pre- and post-treatment RMSE 
#' estimates, and MCMC diagnostics for each of the \eqn{J} models.
leave_donor_out.bscmfit <- function(x, probs = c(0.025, 0.975), ...) {
  stopifnot_(
    checkmate::test_numeric(
      probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values between
     0 and 1."
  )
  stopifnot_(
    !is.null(x$data),
    "The model fit {.arg x} does not contain the original data. You probably 
    used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  outcome <- get_outcome(x)
  treatment <- get_treatment(x)
  donors <- get_donors(x)
  time <- get_time(x)
  times <-get_times(x)
  unit <- get_unit(x)
  data <- x$data
  
  effect <- rmse <- diagnostics <- 
    stats::setNames(vector("list", length(donors)), donors)
  p <- progressr::progressor(along = donors)
  for (donor in donors) {
    p(sprintf(paste0("Estimating the model without donor ", donor, ".")))
    d <- data |> filter(.data[[unit]] != .env$donor)
    fit <- stats::update(
      x, data = d, mcmc_diagnostics = FALSE, save_data = FALSE, ...
    )
    effect[[donor]] <- treatment_effect(fit, probs)
    rmse[[donor]] <- rmse(fit, probs)
    diagnostics[[donor]] <- check_mcmc_diagnostics(fit, warn = FALSE)
  }
  issues <- lapply(diagnostics, \(x) x$has_issues)
  warnifnot_(
    all(!unlist(issues)),
    "Some of the runs resulted in MCMC diagnostic warnings. Check 
    the `diagnostics` element of the output list for details."
  )
  list(effect = effect, rmse = rmse, diagnostics = diagnostics)
}
