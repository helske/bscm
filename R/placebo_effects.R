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
#' @return A list with data frames `effect`, `rmse`, and `diagnostics`, and a
#'   `metadata` list. The data frames contain posterior summaries for each run,
#'   identified by a `placebo` column. For `type = "donor"`, `placebo` is the
#'   treated unit name for the original fit and the donor name for each placebo
#'   run. For `type = "time"`, `placebo` is the assumed treatment start time for
#'   each run. The `metadata` list contains the placebo type, summary
#'   probabilities, model setup, and the placebo labels used. The result can be
#'   visualized with [plot_effects()].
placebo_effects.bscmfit <- function(
  x,
  type,
  L = NULL,
  probs = c(0.025, 0.975),
  ...
) {
  stopifnot_(
    identical(get_N(x), 1L),
    "Placebo effect computation is currently supported only for models with a 
    single treated unit."
  )
  probs <- sort_probs(probs)
  type <- try_(match.arg(type, c("donor", "time")))
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val donor} or {.val time}."
  )
  T_pre <- get_T_pre(x)
  stopifnot_(
    identical(type, "donor") ||
      checkmate::test_integerish(L, len = 1, lower = 2, upper = T_pre - 1),
    "Argument {.arg L} must be a single integer between 2 and {T_pre - 1}, defining 
    the number of time points used for the first fit."
  )
  stopifnot_(
    !is.null(x$data),
    "The model fit {.arg x} does not contain the original data. You probably 
    used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  treatment <- get_treatment(x)
  donors <- get_donors(x)
  time <- get_time(x)
  times <- get_times(x)
  unit <- get_unit(x)
  treated <- get_treated(x)
  data <- x$data
  omega_prior <- get_omega_prior(x)
  if (identical(type, "donor")) {
    placebo_labels <- c(treated, donors)
    effects <- rmses <- diagnostics_list <-
      stats::setNames(vector("list", length(placebo_labels)), placebo_labels)
    effects[[1L]] <- treatment_effect(x, probs = probs)
    rmses[[1L]] <- rmse(x, probs = probs)
    diagnostics_list[[1L]] <- check_mcmc_diagnostics(x, warn = FALSE)
    data <- data |> filter(.data[[unit]] %in% .env$donors)
    end <- times[T_pre]
    p <- progressr::progressor(along = donors)
    for (i in seq_along(donors)) {
      donor <- donors[i]
      p(sprintf(paste0("Estimating the model for donor ", donor, ".")))
      d <- data |>
        mutate(
          "{treatment}" := ifelse(
            .data[[unit]] == .env$donor & .data[[time]] > .env$end,
            1,
            0
          )
        )
      fit <- stats::update(
        x,
        data = d,
        mcmc_diagnostics = FALSE,
        omega_prior = omega_prior,
        ...
      )
      effects[[i + 1L]] <- treatment_effect(fit, probs = probs)
      rmses[[i + 1L]] <- rmse(fit, probs = probs)
      diagnostics_list[[i + 1L]] <- check_mcmc_diagnostics(fit, warn = FALSE)
    }
  }
  if (type == "time") {
    times <- times[L:T_pre]
    placebo_labels <- as.character(times)
    effects <- rmses <- diagnostics_list <-
      stats::setNames(vector("list", length(placebo_labels)), placebo_labels)
    p <- progressr::progressor(along = times)
    for (i in seq_along(times)) {
      p(sprintf(
        paste0("Estimating the model with data up to time ", times[i], ".")
      ))
      d <- data |>
        mutate(
          "{treatment}" := ifelse(
            .data[[unit]] == .env$treated & .data[[time]] > .env$times[i],
            1,
            0
          )
        )
      fit <- stats::update(
        x,
        data = d,
        mcmc_diagnostics = FALSE,
        omega_prior = omega_prior,
        ...
      )
      effects[[i]] <- treatment_effect(fit, probs = probs)
      rmses[[i]] <- rmse(fit, probs = probs)
      diagnostics_list[[i]] <- check_mcmc_diagnostics(fit, warn = FALSE)
    }
  }
  issues <- vapply(diagnostics_list, \(d) d$has_issues, logical(1L))
  warnifnot_(
    all(!issues),
    "Some of the placebo runs resulted in MCMC diagnostic warnings. Check 
    the `diagnostics` element of the output list for details."
  )
  out <- list(
    effect = bind_rows(effects, .id = "placebo"),
    rmse = bind_rows(rmses, .id = "placebo"),
    diagnostics = bind_rows(
      lapply(diagnostics_list, diags2df),
      .id = "placebo"
    ),
    metadata = list(
      type = type,
      probs = probs,
      L = L,
      setup = x$setup,
      placebo = placebo_labels
    )
  )
  class(out) <- "bscm_placebo_effects"
  out
}
