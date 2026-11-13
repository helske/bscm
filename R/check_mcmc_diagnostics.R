#' Check the validity of the posterior of `bscmfit` object
#'
#' This function is automatically called at the end of [bscm()] to check that
#' the output can be trusted in terms of convergence of the MCMC sampling.
#' Checks consists of the common diagnostics of Hamiltonian Monte Carlo variant
#' used by Stan, as well as the Rhat values and effective sample sizes of model
#' parameters. See [rstan::check_hmc_diagnostics()] and
#' [posterior::default_convergence_measures()] for details on the definitions
#' of these.
#'
#' @export
#' @rdname check_mcmc_diagnostics
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param warn \[`logical(1)`]\cr If `TRUE` (the default), generates and
#' prints out a warning in a case of problematic results. Setting
#' this to `FALSE` returns the check results without the warning.
#' @param ... Ignored.
#' @return A list of class `bscmfit_diagnostics` containing the results of the
#' check, see [print.bscmfit_diagnostics()]. If the model has 50 or fewer
#' posterior draws, the diagnostics are not meaningful and are not computed,
#' and `NULL` is returned invisibly with a warning, regardless of `warn`.
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
  
  n_draws <- ndraws(x)
  if (n_draws <= 50L) {
    warning_(
      c(
        "MCMC diagnostics are not meaningful for only few posterior draws.",
        i = "The model contains only {n_draws} draws, so the diagnostics
        were not computed."
      )
    )
    return(invisible(NULL))
  }
  
  n_chains <- nchains(x)
  n_divergences <- rstan::get_num_divergent(get_stanfit(x))
  n_max_treedepth <- rstan::get_num_max_treedepth(get_stanfit(x))
  max_td <- get_stanfit(x)@stan_args[[1L]]$control$max_treedepth %||% 10
  n_low_bfmi <- sum(rstan::get_bfmi(get_stanfit(x)) < 0.2)
  
  sumr <- x |>
    as_draws() |>
    posterior::summarise_draws(posterior::default_convergence_measures())
  idx <- c(
    which.max(sumr$rhat),
    which.min(sumr$ess_bulk),
    which.min(sumr$ess_tail)
  )
  sumr <- dplyr::tibble(
    diagnostic = c("Largest Rhat", "Smallest bulk-ESS", "Smallest tail-ESS"),
    dplyr::as_tibble(sumr[idx, ])
  )
  
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
    x = paste0(
      n_divergences,
      " out of ",
      n_draws,
      " iterations ended with a divergence. This can affect the validity of ",
      "the results. Increasing `adapt_delta` might help."
    ),
    x = paste0(
      "Largest Rhat convergence diagnostic is ",
      signif(sumr$rhat[1], 4),
      " \u2265 1.01 for variable ",
      sumr$variable[1],
      ". As the chains have not necessarily converged, results might not be ",
      "reliable. Increasing the number of iterations might help."
    ),
    `!` = paste0(
      "Smallest effective sample sizes are less than 100 x number of chains ",
      "(bulk-ESS ",
      floor(sumr$ess_bulk[2]),
      " for variable ",
      sumr$variable[2],
      " and tail-ESS ",
      floor(sumr$ess_tail[3]),
      " for ",
      sumr$variable[3],
      "). This affects the accuracy of the posterior summaries. ",
      "Increasing the number of iterations should help."
    ),
    `!` = paste0(
      n_low_bfmi,
      " out of ",
      n_chains,
      " chains had E-BFMI below 0.2, indicating possible issues."
    ),
    i = paste0(
      n_max_treedepth,
      " out of ",
      n_draws,
      " iterations saturated the maximum tree depth of ",
      max_td,
      ". This indicates inefficient sampling, but does not affect the validity",
      " of the results."
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
  out
}
#' Print MCMC Diagnostics
#'
#' Prints whether the diagnostics indicate problems, the messages describing
#' them, and by default the table of the largest Rhat and smallest ESS values.
#'
#' @param x \[`bscmfit_diagnostics`]\cr The diagnostics object returned by
#' [check_mcmc_diagnostics.bscmfit()].
#' @param print_table \[`logical(1)`]\cr If `FALSE`, prints only the overall 
#' result and the messages describing the problems. If `TRUE` (the default), 
#' also prints the table of the largest Rhat and smallest ESS values. 
#' @param ... Ignored.
#' @return input `x` (invisibility).
#' @export
print.bscmfit_diagnostics <- function(x, print_table = TRUE, ...) {
  check_flag(print_table, "print_table")
  if (x$has_issues) {
    cat("\nMCMC diagnostics indicate potential problems:\n")
    cat(cli::format_bullets_raw(x$messages), sep = "\n")
  } else {
    cat("\nMCMC diagnostics indicate no issues.\n")
  }
  if (print_table) {
    cat("\n")
    print_tibble_rows(x$rhat_and_ess)
  }
  invisible(x)
}
