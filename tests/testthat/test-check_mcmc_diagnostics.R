test_that("check_mcmc_diagnostics returns bscmfit_diagnostics", {
  diag <- check_mcmc_diagnostics(fit1_int, warn = FALSE)
  expect_s3_class(diag, "bscmfit_diagnostics")
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
    check_mcmc_diagnostics(fit1_int, warn = "yes"),
    "Argument `warn` must be a single <logical>\\."
  )
})
