#' Print method for bscmfit objects
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @export
print.bscmfit <- function(x, ...) {
  
  setup <- x$setup
  tv <- switch(
    setup$time_varying_effects,
    none = "\n",
    intercept = "time-varying intercept",
    all = "time-varying intercept and regression coefficients"
  )
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Bayesian synthetic control model", deparse(setup$formula))
  if (tv != "\n") {
    cat(" with", tv, "\n")
  }
  cat("Treated unit:", get_treated(x), "\n")
  cat("Number of donors:", length(setup$donors), "\n")
  cat("Number of time periods (pre + post):", setup$T_pre, "+", 
      setup$T_total - setup$T_pre, "\n")
  
  if (!is.null(x$converge)) {
    print(x$converge)
  } else {
    cat(
      paste0(
        "MCMC diagnostics are not available. Use `check_mcmc_diagnostics()` ",
        "on the model fit to obtain these.\n"
      )
    )
  }
  pars <- intersect(
    c("alpha", "beta", "sigma", "sigma_delta"),
    x$stanfit@model_pars
  )
  sumr <- as_draws(x, variable = pars) |> 
    summarise_draws(
      mean, sd, 
      ~ quantile2(.x, probs = c(0.025, 0.975)), 
      default_convergence_measures()
    )
  sumr$variable[sumr$variable == "alpha"] <- "Intercept"
  sumr$variable[sumr$variable == "sigma"] <- "SD_noise"
  sumr$variable[sumr$variable == "sigma_delta"] <- "SD_spline"
  if (length(cf <- x$setup$coef_names)) {
    sumr$variable[grep("^beta\\[", sumr$variable)] <- 
      paste0("Coef_", cf)
    cf <- if (x$setup$has_intercept) cf <- c("Intercept", cf) else cf
    sumr$variable[grep("^sigma_delta\\[", sumr$variable)] <- 
      paste0("SD_spline_", cf)
  }
  cat("\nPosterior summary of main model parameters (excluding weights):\n")
  cat(format(sumr, n = nrow(sumr))[-c(1, 3)], sep = "\n") # remove extra lines
  sumr <- as_draws(x, variable = c("RMSE_pre", "RMSE_post", "R2")) |> 
    summarise_draws(
      mean, sd, 
      ~ quantile2(.x, probs = c(0.025, 0.975)), 
      default_convergence_measures()
    ) |> 
    mutate(
      variable = case_when(
        variable == "RMSE_pre" ~ "Pre-RMSE",
        variable == "RMSE_post" ~ "Post-RMSE",
        variable == "R2" ~ "R^2"
      )
    )
  cat("\nRMSE and Bayesian R^2 values:\n")
  cat(format(sumr, n = nrow(sumr))[-c(1, 3)], sep = "\n") # remove extra lines
  
  invisible(x)
}
