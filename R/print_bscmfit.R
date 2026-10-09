#' Print method for bscmfit objects
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param diagnostics \[`logical(1)`]\cr If `TRUE`, the printed posterior
#'   summaries also contain the MCMC convergence measures and the Monte Carlo
#'   standard error of the mean. Default is `FALSE`, as the convergence of the
#'   sampling is summarised separately, see [check_mcmc_diagnostics()].
#' @param ... Ignored.
#' @rdname print_bscm
#' @export
print.bscmfit <- function(x, diagnostics = FALSE, ...) {
  check_flag(diagnostics, "diagnostics")
  T_pre <- get_T_pre(x)
  T_total <- get_T_total(x)
  N <- get_N(x)
  J <- get_J(x)
  ar1 <- if (x$setup$has_ar1) {
    " with AR(1) residuals"
  } else {
    ""
  }
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Bayesian synthetic control model", deparse(stats::formula(x)), ar1, "\n")
  if (N == 1L) {
    cat("Treated unit:", get_treated(x), "\n")
    cat("Number of donors:", J, "\n")
    cat(
      "Number of time periods (pre + post):",
      T_pre,
      "+",
      T_total - T_pre,
      "\n"
    )
  } else {
    cat("Number of treated units:", N, "\n")
    cat("Number of donors:", J, "\n")
    cat("Number of time periods:", T_total, "\n")
  }
  n_chains <- nchains(x)
  n_warmup <- x$stanfit@sim$warmup
  n_post <- x$stanfit@sim$iter - n_warmup
  cat(
    "MCMC sampling using",
    n_chains,
    "chains with",
    n_warmup,
    "+",
    n_post,
    "iterations took",
    round(max(rowSums(x$elapsed_time$sampling)), 2),
    "seconds\n"
  )
  if (!is.null(x$converge)) {
    print(x$converge)
  } else {
    cat(
      paste0(
        "MCMC diagnostics are not available. Use `check_mcmc_diagnostics()` ",
        "on the model fit object to obtain these.\n"
      )
    )
  }
  print_estimates(x, diagnostics = diagnostics)
  invisible(x)
}

#' Print the posterior summaries shown by `print.bscmfit()`
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#' @param diagnostics \[`logical(1)`]\cr Whether to show the MCMC convergence
#'   measures and the Monte Carlo standard error of the mean.
#' @noRd
print_estimates <- function(x, probs = c(0.025, 0.975), diagnostics = FALSE) {
  single <- get_N(x) == 1L
  pars <- c(
    if (single && has_intercept(x)) "intercept",
    if (has_predictors(x)) "beta",
    if (has_tv_coefs(x) && x$setup$spline_type == "rw2") "nu",
    if (has_tv_coefs(x)) "kappa",
    if (single && has_ar1_error(x)) "rho",
    if (single && has_zeta(x)) "zeta"
  )
  d <- lapply(
    pars,
    \(par) summarise_column(coef_draws(x, par), "coef", probs, diagnostics)
  )
  if (single) {
    d <- c(
      d,
      list(
        sigma(x, probs = probs) |>
          dplyr::mutate(parameter = "Residual SD")
      )
    )
  }
  if (length(d) > 0L) {
    d <- dplyr::bind_rows(d) |>
      dplyr::select(-dplyr::any_of(get_unit(x))) |>
      dplyr::relocate("parameter", .before = 1L) |>
      dplyr::rename(variable = "parameter")
    cat(
      "\nSummary statistics of the time-",
      if (!single) " and unit-",
      "invariant parameters",
      ":\n",
      sep = ""
    )
    print_estimate_table(d, diagnostics)
  }
  if (!single) {
    cat("\nUnit-specific parameters are available via `summary()`.\n")
  }
  if (!is.null(x$data)) {
    treatment <- get_treatment(x)
    e <- list(
      average_treatment_effect(x, average = TRUE, probs = probs) |>
        dplyr::rename(variable = treatment) |>
        dplyr::mutate(
          variable = dplyr::if_else(
            .data$variable == 0,
            "Average pre-treatment effect",
            "Average post-treatment effect"
          )
        )
    )
    if (single) {
      e <- c(
        e,
        list(
          bayes_R2(x, probs = probs) |>
            dplyr::mutate(variable = "Bayesian R2"),
          rmse(x, probs = probs) |>
            dplyr::rename(variable = treatment) |>
            dplyr::mutate(
              variable = dplyr::if_else(
                .data$variable == 0,
                "Pre-treatment RMSE",
                "Post-treatment RMSE"
              )
            )
        )
      )
    }
    e <- dplyr::bind_rows(e) |>
      dplyr::select(-dplyr::any_of(get_unit(x))) |>
      dplyr::relocate("variable", .before = 1L)
    cat(
      "\nAverage treatment effects",
      if (single) " and measures of model fit" else "",
      ":\n",
      sep = ""
    )
    print_estimate_table(e, diagnostics)
  }
}

#' Print a table of posterior summaries without the tibble header lines
#'
#' @param d \[`tibble`]\cr Table of posterior summaries.
#' @param diagnostics \[`logical(1)`]\cr Whether to keep the convergence
#'   measures. The summaries computed by the corresponding methods, such as
#'   [bayes_R2()], always contain them, so they are dropped here.
#' @noRd
print_estimate_table <- function(d, diagnostics) {
  if (!diagnostics) {
    d <- dplyr::select(
      d,
      -dplyr::any_of(c("rhat", "ess_bulk", "ess_tail", "mcse_mean"))
    )
  }
  # remove the dimensions and the column type lines of the tibble output
  cat(format(d, n = nrow(d))[-c(1, 3)], sep = "\n")
}
