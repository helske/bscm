#' @export
#' @rdname leave_donor_out
leave_donor_out <- function(x, ...) {
  UseMethod("leave_donor_out", x)
}
#' Leave-donor-out effects for `bscmfit` 
#' 
#' Re-estimates original model J times by removing one of the J donors from the 
#' model in turn. This can be used to assess how sensitive the results are to the
#' inclusion of specific donors.
#'   
#' @export
#' @rdname leave_donor_out
#' @param x \[`bscmfit`]\cr The output returned by the [bscm()].
#' @param include \[`character()`]\cr Posterior summaries to be computed.   
#' The default is `c("effect", "rmse")`. See [summary.bscmfit()] for details 
#' and list of accepted values. 
#' @param ... Additional parameters passed on to [bscm()].
#' @return A list of data frames containing the estimates.
leave_donor_out.bscmfit <- function(x, include = c("effect", "RMSE"), ...) {
  
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
  
  fits <- stats::setNames(vector("list", length(donors)), donors)
  p <- progressr::progressor(along = donors)
  for (donor in donors) {
    p(sprintf(paste0("Estimating the model without donor ", donor, ".")))
    d <- data |> filter(.data[[unit]] != .env$donor)
    fits[[donor]] <- stats::update(x, data = d, save_data = FALSE, ...)
    fits[[donor]]$convergence <- check_mcmc_diagnostics(
      fits[[donor]], warn = FALSE
    )
    fits[[donor]]$summary <- summary(
      fits[[donor]], include = include
    )
  }
  issues <- names(fits)[unlist(lapply(fits, \(x) x$convergence$has_issues))]
  warnifnot_(
    length(issues) == 0,
    "Some of the placebo runs resulted in MCMC diagnostic warnings. Check the
    diagnostics by `lapply(output, \\(x) x$convergence)`."
  )
  fits
}
