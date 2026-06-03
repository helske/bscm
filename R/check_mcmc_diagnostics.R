#' Check the validity of the posterior of `bscmfit` object
#'
#' This function is automatically called at the end of [bscm()] to check that
#' the output can be trusted in terms of convergence of the MCMC sampling.
#' Checks consists of the common diagnostics of Hamiltonian Monte Carlo variant used by
#' Stan, as well as the Rhat values and effective sample sizes of all model
#' parameters and derived variables. See [rstan::check_hmc_diagnostics()] and
#' [posterior::default_convergence_measures()] for details on the definitions
#' of these.
#'
#' Typical reasons for sampling issues can be a result of short pre-treatment
#' period and/or too many/few donors, units have very different series so that
#' the convex hull assumption is not even approximately plausible. Other
#' reasons include using priors which are incompatible with the data.
#' In addition, too large standard deviation \eqn{kappa} can
#' cause numerical issues (divergences) with Stan.
#'
#' @export
#' @rdname check_mcmc_diagnostics
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param warn \[`logical(1)`]\cr If `TRUE` (the default), generates and
#' (typically) prints out a warning in a case of problematic results. Setting
#' this to `FALSE` silently returns the check results.
#' @param ... Ignored.
#' @return Invisibly returns a list containing the results of the check.
#' @references https://mc-stan.org/learn-stan/diagnostics-warnings.html
check_mcmc_diagnostics <- function(x, ...) {
  UseMethod("check_mcmc_diagnostics", x)
}
#' @export
#' @rdname check_mcmc_diagnostics
#' @examples
#' check_mcmc_diagnostics(fit_single_treated, warn = FALSE)
check_mcmc_diagnostics.bscmfit <- function(x, warn = TRUE, ...) {
  stopifnot_(
    checkmate::test_flag(warn),
    "Argument {.arg warn} must be a single {.cls logical}."
  )
  algorithm <- get_stanfit(x)@stan_args[[1L]]$algorithm
  stopifnot_(
    algorithm %in% c("NUTS", "HMC"),
    "Model was not estimated using MCMC, nothing to diagnose."
  )
  n_chains <- nchains(x)
  n_draws <- ndraws(x)
  stopifnot_(
    n_draws > 50L,
    "MCMC diagnostics are not meaningful for only few posterior draws. 
    The model was estimated with only {n_draws} draws."
  )
  n_divergences <- rstan::get_num_divergent(get_stanfit(x))
  n_max_treedepth <- rstan::get_num_max_treedepth(get_stanfit(x))
  max_td <- get_stanfit(x)@stan_args[[1L]]$control$max_treedepth %||% 10
  n_low_bfmi <- sum(rstan::get_bfmi(get_stanfit(x)) < 0.2)

  sumr <- x |>
    as_draws(parameters = c("y_mean", "y_rep"), include = FALSE) |>
    summarise_draws(posterior::default_convergence_measures())
  idx <- c(
    which.max(sumr$rhat),
    which.min(sumr$ess_bulk),
    which.min(sumr$ess_tail)
  )
  sumr <- cbind(
    data.frame(
      diagnostic = c("Largest Rhat", "Smallest bulk-ESS", "Smallest tail-ESS")
    ),
    sumr[idx, ]
  )
  sumr$rhat <- round(sumr$rhat, 3)
  sumr$ess_bulk <- round(sumr$ess_bulk)
  sumr$ess_tail <- round(sumr$ess_tail)
  rhat_ok <- sumr$rhat[1] < 1.01

  ess_bulk_ok <- sumr$ess_bulk[2] > (n_chains * 100)
  ess_tail_ok <- sumr$ess_tail[3] > (n_chains * 100)

  issues <- c(
    divergences = n_divergences > 0,
    rhat = !rhat_ok,
    ess = !ess_bulk_ok || !ess_tail_ok,
    n_low_bfmi = n_low_bfmi > 0,
    max_treedepth = n_max_treedepth > 0
  )
  messages <- c(
    x = paste(
      n_divergences,
      "out of",
      n_draws,
      "iterations ended with a divergence. This can affect the validity of",
      "the results. Increasing `adapt_delta` might help."
    ),
    x = paste0(
      "Largest Rhat convergence diagnostic is ",
      sumr$rhat[1],
      " > 1.01 for variable ",
      sumr$variable[1],
      ". As the chains have not necessarily converged, results might not be ",
      "reliable. Increasing the number of iterations might help."
    ),
    `!` = paste0(
      "Smallest effective samples sizes are less than 100 x number of chains ",
      "(bulk-ESS ",
      sumr$ess_bulk[2],
      " for variable ",
      sumr$variable[2],
      " and tail-ESS ",
      sumr$ess_tail[3],
      " for ",
      sumr$variable[3],
      "). This affects the accuracy of the posterior summaries. ",
      "Increasing the number of iterations might help."
    ),
    `!` = paste(
      n_low_bfmi,
      "out of",
      n_chains,
      "chains had E-BFMI below 0.2, indicating possible issues."
    ),
    i = paste(
      n_max_treedepth,
      "out of",
      n_draws,
      "iterations saturated the maximum tree depth of",
      max_td,
      "This indicates inefficient sampling, but does not affect the validity",
      "of the results."
    )
  )[which(issues)]

  out <- dplyr::lst(
    n_divergences,
    n_max_treedepth,
    n_low_bfmi,
    rhat_and_ess = sumr,
    has_issues = any(issues),
    messages = messages
  )
  if (warn && any(issues)) {
    cli::cli_warn(messages)
  }
  class(out) <- "bscmfit_diagnostics"
  invisible(out)
}
#' Print MCMC Diagnostics
#'
#' @param x \[`bscmfit_diagnostics`]\cr The diagnostics object returned by
#' [check_mcmc_diagnostics.bscmfit()].
#' @param print_table \[`logical(1)`] If `NULL` (the default) only prints the
#' table of largest Rhat and smallest ESS values in case diagnostics indicate
#' problems. If `TRUE` or `FALSE`, always prints or does not print the table.
#' @param ... Ignored.
#' @return input `x` (invisibility).
#' @export
print.bscmfit_diagnostics <- function(x, print_table = NULL, ...) {
  stopifnot_(
    checkmate::test_flag(print_table, null.ok = TRUE),
    "Argument {.arg print_table} must be a single {.cls logical} value."
  )
  if (x$has_issues) {
    if (is.null(print_table)) {
      print_table <- TRUE
    }
    warning_(c("\nMCMC diagnostics indicate potential problems:", x$messages))
  } else {
    if (is.null(print_table)) {
      print_table <- FALSE
    }
    cat("\nMCMC diagnostics indicate no issues. \n")
  }
  if (print_table) {
    print(x$rhat_and_ess)
  }
  invisible(x)
}
