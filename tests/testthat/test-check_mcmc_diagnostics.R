test_that("check_mcmc_diagnostics returns bscmfit_diagnostics", {
  suppressWarnings(diag <- check_mcmc_diagnostics(fit_intercept, warn = FALSE))
  expect_s3_class(diag, "bscmfit_diagnostics")
  expect_gt(length(diag$me), 0)
  expected <- c(
    "n_divergences", "n_max_treedepth", "n_low_bfmi",
    "rhat_and_ess", "has_issues", "messages"
  )
  expect_true(all(expected %in% names(diag)))
  expect_equal(
    diag$rhat_and_ess$diagnostic,
    c("Largest Rhat", "Smallest bulk-ESS", "Smallest tail-ESS")
  )
  expect_error(
    check_mcmc_diagnostics(fit_intercept, check_all = "yes"),
    "Argument `check_all` must be a single <logical>\\."
  )
  expect_error(
    check_mcmc_diagnostics(fit_intercept, warn = "yes"),
    "Argument `warn` must be a single <logical>\\."
  )
})
