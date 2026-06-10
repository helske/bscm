#' Print method for bscmfit objects
#'
#' @rdname print_bscm
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return Returns `x` (invisibly).
#' @export
print.bscmfit <- function(x, ...) {
  T_pre <- get_T_pre(x)
  T_total <- get_T_total(x)
  N <- get_N(x)
  J <- get_J(x)
  ar1 <- if (x$setup$error == "ar1") {
    " with AR(1) residuals"
  } else ""
  
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
    "chains, each with",
    n_warmup,
    "+",
    n_post,
    "iterations took",
    round(max(rowSums(x$elapsed_time$sampling)), 2),
    "seconds for the slowest chain\n"
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
  invisible(x)
}

#' @rdname print_bscm
#' @param x \[`bscmfit` or `summary_bscmfit`]\cr Output from [bscm()] or
#' [summary.bscmfit()].
#' @param ... Ignored.
#' @return Returns `x` (invisibly).
#' @export
print.summary_bscmfit <- function(x, ...) {
  cat("\n")
  cat(format(x, n = nrow(x))[-c(1, 3)], sep = "\n") # remove extra lines
  invisible(x)
}
