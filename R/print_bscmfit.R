#' Print method for bscmfit objects
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#'   Default is `c(0.025, 0.975)`.
#' @param diagnostics \[`logical(1)`]\cr If `TRUE` (the default), 
#'   the printed posterior summaries contain the MCMC convergence measures and 
#'   the Monte Carlo standard error of the mean.
#' @param ... Ignored.
#' @rdname print_bscm
#' @export
print.bscmfit <- function(x, probs = c(0.025, 0.975), diagnostics = TRUE, ...) {
  probs <- sort_probs(probs)
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
  print_estimates(x, probs = probs, diagnostics = diagnostics)
  invisible(x)
}

#' Print the posterior summaries shown by `print.bscmfit()`
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries.
#' @param diagnostics \[`logical(1)`]\cr Whether to show the MCMC convergence
#'   measures and the Monte Carlo standard error of the mean.
#' @noRd
print_estimates <- function(x, probs, diagnostics) {
  single <- get_N(x) == 1L
  d <- summary(x, probs = probs, diagnostics = diagnostics)
  if (!single) {
    d <- dplyr::filter(d, .data$parameter %in% c("coef", "tv_sd", "tv_slope"))
  }
  par <- d$parameter
  term <- d[["term"]] %||% rep(NA_character_, nrow(d))
  d$variable <- dplyr::case_when(
    par == "coef" ~ term,
    par == "tv_sd" ~ paste0("sd(", term, ")"),
    par == "tv_slope" ~ paste0("slope(", term, ")"),
    .default = par
  )
  sections <- list(
    "Coefficients" = c("intercept", "coef"),
    "Time variation of coefficients" = c("tv_sd", "tv_slope"),
    "Error" = c("sigma", "ar1")
  )
  for (section in names(sections)) {
    rows <- d[d$parameter %in% sections[[section]], ]
    if (nrow(rows) > 0L) {
      cat("\n", section, ":\n", sep = "")
      rows |>
        dplyr::select(-dplyr::any_of(c("parameter", "term", get_unit(x)))) |>
        dplyr::relocate("variable", .before = 1L) |>
        print_estimate_table()
    }
  }
  if (!single) {
    cat("\nUnit-specific parameters are available via `summary()`.\n")
  }
  cat("\nDonor weights:\n")
  effective_donors(
    x,
    average = TRUE,
    probs = probs,
    diagnostics = diagnostics
  ) |>
    dplyr::mutate(
      variable = if (single) {
        "Effective number of donors"
      } else {
        "Average effective number of donors"
      }
    ) |>
    dplyr::select(-dplyr::any_of(get_unit(x))) |>
    dplyr::relocate("variable", .before = 1L) |>
    print_estimate_table()
  if (!is.null(x$data)) {
    treatment <- get_treatment(x)
    e <- list(
      average_treatment_effect(
        x, 
        average = TRUE, 
        probs = probs, 
        diagnostics = diagnostics
      ) |>
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
          bayes_R2(x, probs = probs, diagnostics = diagnostics) |>
            dplyr::mutate(variable = "Bayesian R2"),
          rmse(x, probs = probs, diagnostics = diagnostics) |>
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
    print_estimate_table(e)
  } else {
    cat(
      "\nSome of the summary statistics are not available because the ",
      "original data are not stored in the model fit object.\n",
      sep = ""
    )
  }
}

#' Print a table of posterior summaries without the tibble header lines
#'
#' @param d \[`tibble`]\cr Table of posterior summaries.
#' @noRd
print_estimate_table <- function(d) {
  # remove the dimensions and the column type lines of the tibble output
  cat(format(d, n = nrow(d), width = Inf)[-c(1, 3)], sep = "\n")
}
