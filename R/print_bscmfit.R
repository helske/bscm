#' Print method for bscmfit objects
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return The summary output of the model.
#' @export
print.bscmfit <- function(x, ...) {
  
  setup <- x$setup
  T_pre <- setup$T_pre
  T_total <- setup$T_total
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Bayesian synthetic control model", deparse(stats::formula(x)))
  if (setup$time_varying_effects) {
    cat(" with time-varying coefficients \n")
  } else {
    cat("\n")
  }
  cat("Treated unit:", get_treated(x), "\n")
  cat("Number of donors:", length(setup$donors), "\n")
  cat("Number of time periods (pre + post):", T_pre, "+", T_total - T_pre, "\n")
  
  cat(
    "MCMC sampling time:", max(rowSums(x$elapsed_time$sampling)), "seconds\n"
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
  pars <- c(
    if (setup$has_intercept) "alpha", if (length(setup$predictors) > 0L) "beta", 
    if (setup$time_varying_effects) "tau",
    "sigma", paste0("effect[", T_pre + 1, "]"), "avg_effect_post",
    "effective_donors", "RMSE_pre", "RMSE_post", "R2"
  )
  timevar <- setup$time
  t0 <- setup$times[T_pre + 1]
  beta_map <- stats::setNames(
    paste0("Coef_", x$setup$coef_names),
    paste0("beta[", seq_along(x$setup$coef_names), "]")
  )
  
  sumr <- as_draws(x, parameters = pars) |> 
    summarise_draws(
      mean, sd, 
      ~ quantile2(.x, probs = c(0.025, 0.975)), 
      default_convergence_measures()
    ) |> 
    mutate(
      variable = case_when(
        variable == "alpha" ~ "Intercept",
        variable %in% names(beta_map) ~ beta_map[variable],
        variable == "tau" ~ "Base SD of the time-varying imbalance effects",
        variable == "sigma" ~ "Residual SD",
        variable == paste0("effect[", T_pre + 1, "]") ~ 
          paste("Treatment effect at", timevar, t0),
        variable == "avg_effect_post" ~ "Average treatment effect",
        variable == "effective_donors" ~ "Number of effective donors",
        variable == "RMSE_pre" ~ "Pre-RMSE",
        variable == "RMSE_post" ~ "Post-RMSE",
        variable == "R2" ~ "Bayesian R^2"
      )
    )
  cat("\n")
  cat(format(sumr, n = nrow(sumr))[-c(1, 3)], sep = "\n") # remove extra lines
  invisible(sumr)
}
