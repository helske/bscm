#' @export
#' @rdname placebo_effects
placebo_effects <- function(x, ...) {
  UseMethod("placebo_effects", x)
}
#' Placebo effects for `bscmfit` 
#' 
#' For the in-space placebo, original model is re-estimated using each donor as 
#' the treated unit in turn (ignoring the original, true treated). For the 
#' in-time placebo, we still estimate the treatment effect for the original 
#' treated, but consider only the pretreatment data, and moving the start of the 
#' treatment from \eqn{L + 1} to \eqn{T_0 - 1}, where \eqn{L} is the minimum 
#' number of pre-treatment time points to use, and \eqn{T_0} is the time of the 
#' true intervention.
#'   
#' @export
#' @rdname placebo_effects
#' @param x \[`bscmfit`]\cr The output returned by the [bscm()].
#' @param type \[`character(1)`]\cr Type of the placebo effects to compute. 
#' Either `"space"` for in-space placebos or `"time"` for in-time placebos. See 
#' details.
#' @param L \[`integer(1)`]\cr If `type = "time`, minimum number of observations 
#' to use for the in-time placebos, i.e. the number of pre-treatment time points 
#' for the first fit. For too small `L`, estimation can be unstable, so you 
#' should likely use at least `L = 10`.
#' @param summaries \[`character()`]\cr Posterior summaries to be computed.   
#' The default is `c("effect", "rmse")`. See [summary.bscmfit()] for details 
#' and list of accepted values. 
#' @param ... Additional parameters passed on to [bscm()].
#' @return A list of data frames containing the estimates.
placebo_effects.bscmfit <- function(x, type, L = NULL, 
                                    summaries = c("effect", "RMSE"), ...) {
  
  type <- try_(match.arg(type, c("space", "time")))
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val space} or {.val time}."
  )
  T_pre <- get_T_pre(x)
  stopifnot_(
    identical(type, "space") || 
      checkmate::test_integerish(L, len = 1, lower = 2, upper = T_pre - 1),
    "Argument {.arg L} must be a single integer between 2 and {T_pre}, defining 
    the number of time points for the first fit."
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
  if (identical(type, "space")) {
    fits <- stats::setNames(vector("list", length(donors)), donors)
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
      fits[[donor]] <- stats::update(
        x, data = d, save_data = FALSE, ...
      )
    }
  }
  if (type == "time") {
    end <- times[T_pre]
    treated <- x$setup$treated
    data <- data |> 
      filter(.data[[time]] <= .env$end)
    times <- times[(L + 1):T_pre]
    fits <- stats::setNames(vector("list", T_pre - L), times)
    p <- progressr::progressor(along = fits)
    for (t in times) {
      p(sprintf(paste0("Estimating the model for up to time ", t, ".")))
      d <- data |> 
        filter(.data[[time]] <= .env$t) |> 
        mutate(
          "{treatment}" := ifelse(
            .data[[unit]] == .env$treated & .data[[time]] == .env$t, 1, 0
          )
        )
      fits[[as.character(t)]] <- stats::update(
        x, data = d, save_data = FALSE, ...
      )
    }
  }
  issues <- names(fits)[unlist(lapply(fits, \(x) x$convergence$has_issues))]
  warnifnot_(
    length(issues) == 0,
    "Some of the placebo runs resulted in MCMC diagnostic warnings. Check the
    diagnostics by `lapply(output, \\(x) x$convergence)`."
  )
  fits
}
