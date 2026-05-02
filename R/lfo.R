#' Leave-Future-Out Cross-Validation
#'
#' @export
#' @rdname lfo
lfo <- function(x, ...) {
  UseMethod("lfo", x)
}

#' Leave-future-out (LFO) cross-validation for Bayesian synthetic control models
#'
#' Estimates the leave-future-out (LFO) expected log predictive density (ELPD)
#' for `bscmfit` models. The LFO-CV is performed over the pre-treatment period:
#' at each step `t` (from `L` to `T_pre_min - 1`), the model is evaluated on
#' its ability to predict the next pre-treatment observation `y[t+1]` having
#' been fitted on `y[1:t]`. For models with multiple treated units, all units
#' must have at least `L + 1` pre-treatment observations, and only the first
#' `min(T_pre)` time points are used so the number of treated units stays
#' constant throughout the evaluation period. In this case, a joint ELPD of
#' all treated units is computed at each step.
#'
#' Two methods are available:
#'
#' * **Exact LFO** (`exact = TRUE`): the model is refitted at every step.
#'   This is exact but expensive.
#'
#' * **PSIS-LFO** (`exact = FALSE`, the default): uses Pareto-smoothed
#'   importance sampling (PSIS) to avoid refitting at every step. The model is
#'   only refitted when the Pareto `k` diagnostic exceeds `k_thres`. At refit
#'   points the ELPD is computed exactly; at other points it is approximated
#'   by PSIS.
#'
#' @references Paul-Christian Bürkner, Jonah Gabry & Aki Vehtari (2020).
#'   Approximate leave-future-out cross-validation for Bayesian time series
#'   models. *Journal of Statistical Computation and Simulation*, 90(14),
#'   2499--2523. <https://doi.org/10.1080/00949655.2020.1783262>
#'
#' @export
#' @rdname lfo
#' @param x \[`bscmfit`]\cr The output returned by [bscm()].
#' @param L \[`integer(1)`]\cr Minimum number of pre-treatment observations
#'   used for the first fit. Must satisfy `2 <= L <= T_pre_min - 2`, where
#'   `T_pre_min = min(T_pre)`. Too small value of `L` can lead to unstable
#'   estimation, so a value of at least `10` is recommended.
#' @param exact \[`logical(1)`]\cr If `TRUE`, computes exact LFO by refitting
#'   at every step. If `FALSE` (the default), uses the approximate PSIS-LFO
#'   method.
#' @param k_thres \[`numeric(1)`]\cr Threshold for the Pareto `k` diagnostic
#'   that triggers a model refit. Default is `0.7`. Ignored when `exact = TRUE`.
#' @param ... Additional arguments passed on to [bscm()] when refitting.
#' @return An object of class `bscm_lfo`, a list with components:
#'   * `ELPD`: Total expected log predictive density.
#'   * `ELPDs`: Vector of per-step ELPDs (length `T_pre_min - L`). Element `k`
#'     is the ELPD for predicting observation `L + k`.
#'   * `ks`: Pareto `k` values (length `T_pre_min - L - 1`, `NULL` for exact
#'     LFO).
#'   * `refits`: Time indices at which the model was re-estimated.
#'   * `L`: The value of `L` used.
#'   * `T_pre_min`: The minimum pre-treatment period length across treated
#'     units.
#'   * `times`: Vector of all unique time values from the original data.
#'   * `time_var`: Name of the time variable.
#'   * `k_thres`: The Pareto k threshold used.
#' @examples
#' \dontrun{
#' fit <- bscm(
#'   y ~ 1, single_treated, "treatment", "time", "id",
#'   chains = 1, cores = 1, refresh = 0
#' )
#' lfo_approx <- lfo(fit, L = 10)
#' lfo_approx
#'
#' lfo_exact <- lfo(fit, L = 10, exact = TRUE)
#' lfo_exact
#' }
lfo.bscmfit <- function(
  x,
  L,
  exact = FALSE,
  k_thres = 0.7,
  ...
) {
  stopifnot_(
    !is.null(x$data),
    "LFO requires the original data. Refit the model with
    {.code save_data = TRUE}."
  )
  T_pre <- get_T_pre(x)
  T_pre_min <- min(T_pre)
  stopifnot_(
    checkmate::test_integerish(L, len = 1L, lower = 2L, upper = T_pre_min - 2L),
    "Argument {.arg L} must be a single integer between 2 and {T_pre_min - 2}."
  )

  time <- get_time(x)
  times <- get_times(x)
  unit <- get_unit(x)
  treatment <- get_treatment(x)
  treated <- get_treated(x)
  omega_prior <- get_omega_prior(x)

  refit_at <- function(t_lfo) {
    d <- x$data |>
      mutate(
        "{treatment}" := ifelse(
          .data[[unit]] %in% .env$treated & .data[[time]] > .env$times[t_lfo],
          1L,
          0L
        )
      )
    stats::update(
      x,
      data = d,
      compute_predictions = FALSE,
      mcmc_diagnostics = FALSE,
      save_data = FALSE,
      omega_prior = omega_prior,
      refresh = 0,
      ...
    )
  }

  treated_data <- x$data |> filter(.data[[unit]] %in% .env$treated)
  n_steps <- T_pre_min - L
  elpds <- rep(NA_real_, n_steps)
  p <- progressr::progressor(steps = n_steps)

  if (exact) {
    refits <- L:(T_pre_min - 1L)
    ks <- NULL
    for (t in L:(T_pre_min - 1L)) {
      p(sprintf("Fitting model with %d pre-treatment observations.", t))
      fit_t <- refit_at(t)
      ll_next <- log_lik_at_times(fit_t, treated_data, t + 1L)
      elpds[t - L + 1L] <- log_mean_exp(ll_next[, 1L])
    }
  } else {
    p(sprintf("Fitting model with %d pre-treatment observations.", L))
    fit_past <- refit_at(L)
    ll_from_refit <- log_lik_at_times(
      fit_past,
      treated_data,
      seq(L + 1L, T_pre_min)
    )

    elpds[1L] <- log_mean_exp(ll_from_refit[, 1L])
    i_refit <- L
    refits <- L
    ks <- rep(NA_real_, n_steps - 1L)

    for (step in seq(2L, n_steps)) {
      t <- L + step - 1L
      n_bridge <- t - i_refit
      logratio <- rowSums(ll_from_refit[, seq_len(n_bridge), drop = FALSE])
      psis_obj <- suppressWarnings(loo::psis(logratio))
      k <- loo::pareto_k_values(psis_obj)
      ks[step - 1L] <- k
      if (k > k_thres) {
        i_refit <- t
        refits <- c(refits, t)
        p(sprintf("Refitting model with %d pre-treatment observations.", t))
        fit_past <- refit_at(t)
        ll_from_refit <- log_lik_at_times(
          fit_past,
          treated_data,
          seq(t + 1L, T_pre_min)
        )
        elpds[step] <- log_mean_exp(ll_from_refit[, 1L])
      } else {
        p(sprintf("PSIS approximation at step %d (k = %.2f).", t, k))
        lw <- loo::weights.importance_sampling(psis_obj, normalize = TRUE)[, 1L]
        elpds[step] <- log_sum_exp(lw + ll_from_refit[, step - (i_refit - L)])
      }
    }
  }

  out <- list(
    ELPD = sum(elpds),
    ELPD_SE = stats::sd(elpds) * sqrt(length(elpds)),
    ELPDs = elpds,
    ks = ks,
    refits = refits,
    L = L,
    T_pre_min = T_pre_min,
    times = times,
    time_var = time,
    k_thres = k_thres
  )
  class(out) <- "bscm_lfo"
  out
}
#' Print method for LFO cross-validation output
#'
#' Prints the summary of the leave-future-out cross-validation.
#'
#' @param x \[`bscm_lfo`]\cr Output from [lfo()].
#' @param ... Ignored.
#' @return Returns `x` invisibly.
#' @export
print.bscm_lfo <- function(x, ...) {
  approx <- if (!is.null(x$ks)) "Approximate" else "Exact"
  cat(
    "\n",
    approx,
    " LFO starting from ",
    x$time_var,
    " = ",
    x$times[x$L],
    sep = ""
  )
  if (!is.null(x$ks)) {
    cat(
      "\nModel was re-estimated at ",
      x$time_var,
      " = ",
      paste(x$times[x$refits], collapse = ", "),
      " (Based on Pareto k threshold of ",
      x$k_thres,
      ")\n",
      sep = ""
    )
  } else {
    cat("\nModel was re-estimated at every time step (exact LFO)\n")
  }
  cat("\nEstimated expected log predictive density (ELPD):", x$ELPD)
  cat("\nStandard error estimate of the ELPD:", x$ELPD_SE)
  invisible(x)
}

#' Plot Pareto k diagnostics from LFO cross-validation
#'
#' Plots Pareto k values against the corresponding time variable values,
#' with a horizontal line at the refitting threshold. Only available (and relevant)
#' for PSIS-LFO (not exact LFO).
#'
#' @param x \[`bscm_lfo`]\cr Output from [lfo()].
#' @param ... Ignored.
#' @return A `ggplot` object, or `NULL` invisibly for exact LFO.
#' @export
plot.bscm_lfo <- function(x, ...) {
  if (is.null(x$ks)) {
    message(
      "No Pareto k values available (output is based on exact LFO). Nothing to plot."
    )
    return(invisible(NULL))
  }
  d <- data.frame(
    k = x$ks,
    time = x$times[x$L + seq_len(length(x$ks))]
  )
  d$threshold <- d$k > x$k_thres
  # avoid NSE notes from R CMD check
  time <- k <- threshold <- NULL
  ggplot(d, aes(x = time, y = k)) +
    geom_point(
      aes(color = threshold),
      shape = 3,
      show.legend = FALSE,
      alpha = 0.5
    ) +
    geom_hline(
      yintercept = x$k_thres,
      linetype = 2,
      color = "red2"
    ) +
    scale_color_manual(values = c("cornflowerblue", "darkblue")) +
    labs(x = x$time_var, y = "Pareto k")
}
