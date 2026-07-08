test_that("check_mcmc_diagnostics returns bscmfit_diagnostics", {
  diag <- check_mcmc_diagnostics(fit1_int, warn = FALSE)
  expect_s3_class(diag, "bscmfit_diagnostics")
  expected <- c(
    "n_divergences",
    "n_max_treedepth",
    "n_low_bfmi",
    "rhat_and_ess",
    "has_issues",
    "messages"
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

test_that("print.bscmfit_diagnostics returns x invisibly", {
  diag <- check_mcmc_diagnostics(fit1_int, warn = FALSE)
  expect_invisible(print(diag))
})

test_that("print.bscmfit_diagnostics indicates no issues for a well-fit model", {
  diag <- check_mcmc_diagnostics(fit1_int, warn = FALSE)
  out <- capture.output(print(diag))
  expect_true(any(grepl("no issues", out)))
})

test_that("print.bscmfit_diagnostics with print_table = TRUE includes Rhat row", {
  diag <- check_mcmc_diagnostics(fit1_int, warn = FALSE)
  out <- capture.output(print(diag, print_table = TRUE))
  expect_true(any(grepl("Rhat", out, fixed = TRUE)))
})

test_that("print.bscmfit_diagnostics validates print_table argument", {
  diag <- check_mcmc_diagnostics(fit1_int, warn = FALSE)
  expect_error(print(diag, print_table = "yes"), "single <logical>")
})
