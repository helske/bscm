#' @export
#' @rdname placebo_effects
placebo_effects <- function(x, ...) {
  UseMethod("placebo_effects", x)
}
#' Placebo effects of a Bayesian synthetic control model
#' 
#' For the in-space placebo (`type = "donor"`), original model is re-estimated 
#' using each donor as the treated unit in turn (omitting the original, true 
#' treated unit). The obtained effect estimates should be around zero, under 
#' the assumption that the treatment affected only the true treated unit.
#' 
#' For the in-time placebo (`type = "time"`) we still 
#' estimate the treatment effect for the original treated, but move the start 
#' of the treatment from \eqn{L + 1} to \eqn{T_{pre} + 1}, where \eqn{L} is the 
#' minimum number of pre-treatment time points to use, and \eqn{T_{pre}} is the 
#' true last pre-treatment time point. In all cases, the obtained treatment 
#' effects should fluctuate around zero for time points before the true 
#' treatment time, under the assumption of no anticipation effects. 
#'   
#' @export
#' @rdname placebo_effects
#' @param x \[`bscmfit`]\cr The output returned by the [bscm()].
#' @param type \[`character(1)`]\cr Type of the placebo effects to compute. 
#' Either `"donor"` for in-space placebos, `"time"` for in-time placebos. See 
#' details.
#' @param L \[`integer(1)`]\cr If `type = "time`, minimum number of observations 
#' to use for the in-time placebos, i.e. the number of pre-treatment time points 
#' for the first fit. For too small `L`, estimation can be unstable, so you 
#' should likely use at least `L = 10` or so.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries of the 
#' treatment effects and RMSE estimates. Default is `c(0.025, 0.975)`.
#' @param ... Additional arguments passed on to [bscm()].
#' @return A list with three elements: `effect`, `rmse`, and `diagnostics`, 
#' containing the treatment effect estimates, pre- and post-treatment RMSE 
#' estimates, and MCMC diagnostics for each placebo run.
placebo_effects.bscmfit <- function(x, type, L = NULL, 
                                    probs = c(0.025, 0.975), ...) {
  stopifnot_(
    identical(get_N(x), 1L),
    "Placebo effect computation is currently supported only for models with a 
    single treated unit."
  )
  test_probs(probs)
  type <- try_(match.arg(type, c("donor", "time")))
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val donor} or {.val time}."
  )
  T_pre <- get_T_pre(x)
  stopifnot_(
    identical(type, "donor") || 
      checkmate::test_integerish(L, len = 1, lower = 2, upper = T_pre - 1),
    "Argument {.arg L} must be a single integer between 2 and {T_pre}, defining 
    the number of time points used for the first fit."
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
  times <- get_times(x)
  unit <- get_unit(x)
  data <- x$data
  if (identical(type, "donor")) {
    effect <- rmse <- diagnostics <- 
      stats::setNames(vector("list", length(donors)), donors)
    data <- data |> filter(.data[[unit]] %in% .env$donors)
    end <- times[T_pre]
    p <- progressr::progressor(along = donors)
    for (donor in donors) {
      p(sprintf(paste0("Estimating the model for donor ", donor, ".")))
      d <- data |> 
        mutate(
          "{treatment}" := ifelse(
            .data[[unit]] == .env$donor & .data[[time]] > .env$end, 1, 0
          )
        )
      fit <- stats::update(
        x, data = d, mcmc_diagnostics = FALSE, save_data = FALSE, ...
      )
      effect[[donor]] <- treatment_effect(fit, probs)
      rmse[[donor]] <- rmse(fit, probs)
      diagnostics[[donor]] <- check_mcmc_diagnostics(fit, warn = FALSE)
    }
  }
  if (type == "time") {
    treated <- x$setup$treated
    times <- times[L:T_pre]
    effect <- rmse <- diagnostics <- 
      stats::setNames(vector("list", T_pre - L + 1), times)
    p <- progressr::progressor(along = times)
    for (t in times) {
      p(sprintf(
        paste0("Estimating the model with data up to time ", t, ".")
      ))
      d <- data |> 
        mutate(
          "{treatment}" := ifelse(
            .data[[unit]] == .env$treated & .data[[time]] > .env$t, 1, 0
          )
        )
      fit <- stats::update(
        x, data = d, mcmc_diagnostics = FALSE, save_data = FALSE, ...
      )
      idx <- as.character(t)
      effect[[idx]] <- treatment_effect(fit, probs)
      rmse[[idx]] <- rmse(fit, probs)
      diagnostics[[idx]] <- check_mcmc_diagnostics(fit, warn = FALSE)
    }
  }
  issues <- lapply(diagnostics, \(x) x$has_issues)
  warnifnot_(
    all(!unlist(issues)),
    "Some of the placebo runs resulted in MCMC diagnostic warnings. Check 
    the `diagnostics` element of the output list for details."
  )
  list(effect = effect, rmse = rmse, diagnostics = diagnostics)
}
