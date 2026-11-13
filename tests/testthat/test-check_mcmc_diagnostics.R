test_that("check_mcmc_diagnostics returns bscmfit_diagnostics", {
  diag <- check_mcmc_diagnostics(fit_single_treated, warn = FALSE)
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
  expect_s3_class(diag$rhat_and_ess, "tbl_df")
  expect_equal(
    diag$rhat_and_ess$diagnostic,
    c("Largest Rhat", "Smallest bulk-ESS", "Smallest tail-ESS")
  )
  expect_error(
    check_mcmc_diagnostics(fit1_int, warn = "yes"),
    "Argument `warn` must be a single <logical>\\."
  )
})

test_that("check_mcmc_diagnostics returns the results visibly", {
  expect_visible(check_mcmc_diagnostics(fit_single_treated, warn = FALSE))
})

test_that("print.bscmfit_diagnostics returns x invisibly", {
  diag <- check_mcmc_diagnostics(fit_single_treated, warn = FALSE)
  capture.output(expect_invisible(print(diag)))
})

test_that("print.bscmfit_diagnostics indicates no issues for a well-fit model", {
  diag <- check_mcmc_diagnostics(fit_single_treated, warn = FALSE)
  out <- capture.output(print(diag))
  expect_true(any(grepl("no issues", out)))
})

test_that("print.bscmfit_diagnostics shows the table by default", {
  diag <- check_mcmc_diagnostics(fit_single_treated, warn = FALSE)
  out <- capture.output(print(diag))
  expect_true(any(grepl("Largest Rhat", out, fixed = TRUE)))
  out <- capture.output(print(diag, print_table = FALSE))
  expect_false(any(grepl("Largest Rhat", out, fixed = TRUE)))
})

test_that("print.bscmfit_diagnostics prints the problems without warning", {
  diag <- check_mcmc_diagnostics(fit_single_treated, warn = FALSE)
  diag$has_issues <- TRUE
  diag$messages <- c(x = "First problem.", i = "Second problem.")
  expect_no_warning(out <- capture.output(print(diag)))
  expect_true(any(grepl("potential problems", out, fixed = TRUE)))
  expect_true(any(grepl("First problem.", out, fixed = TRUE)))
  expect_true(any(grepl("Second problem.", out, fixed = TRUE)))
})

test_that("print.bscmfit_diagnostics validates print_table argument", {
  diag <- check_mcmc_diagnostics(fit_single_treated, warn = FALSE)
  expect_error(print(diag, print_table = "yes"), "single <logical>")
})

test_that("check_mcmc_diagnostics warns and returns NULL for few draws", {
  fit <- suppressWarnings(
    bscm(
      y ~ 1,
      data = single_treated,
      chains = 1,
      iter = 10,
      refresh = 0,
      seed = 1,
      mcmc_diagnostics = FALSE
    )
  )
  expect_warning(
    d <- check_mcmc_diagnostics(fit),
    "not meaningful for only few posterior draws"
  )
  expect_null(d)
  expect_warning(
    check_mcmc_diagnostics(fit, warn = FALSE),
    "not meaningful for only few posterior draws"
  )
})
