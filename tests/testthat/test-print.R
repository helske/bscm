test_that("print.bscmfit returns x invisibly", {
  out <- capture.output(result <- print(fit_single_treated))
  expect_identical(result, fit_single_treated)
})

test_that("print.bscmfit output contains expected header", {
  out <- capture.output(print(fit_single_treated))
  expect_true(any(grepl("Bayesian synthetic control model", out)))
})

test_that("print.bscmfit output contains MCMC info", {
  out <- capture.output(print(fit_single_treated))
  expect_true(any(grepl("MCMC sampling", out)))
})

test_that("print.bscmfit shows the diagnostics details only with issues", {
  out <- capture.output(print(fit_single_treated))
  expect_true(any(grepl("no issues", out, fixed = TRUE)))
  expect_false(any(grepl("Largest Rhat", out, fixed = TRUE)))
})

test_that("print.bscmfit shows treated unit info for single treated", {
  out <- capture.output(print(fit_single_treated))
  expect_true(any(grepl("Treated unit:", out)))
  expect_true(any(grepl("Number of donors:", out)))
})

test_that("print.bscmfit shows number of treated units for N > 1", {
  suppressWarnings(out <- capture.output(print(fitN_tv)))
  expect_true(any(grepl("Number of treated units:", out)))
})

test_that("print.bscmfit shows AR(1) residuals for AR(1) model", {
  suppressWarnings(out <- capture.output(print(fit1_noint)))
  expect_true(any(grepl("AR(1) residuals", out, fixed = TRUE)))
})

test_that("print.bscmfit does not show AR(1) for IID model", {
  out <- capture.output(print(fit_single_treated))
  expect_false(any(grepl("AR(1) residuals", out, fixed = TRUE)))
})

test_that("print.bscmfit shows the summaries of the model parameters", {
  suppressWarnings(out <- capture.output(print(fit1_x)))
  expect_true(all(c("Coefficients:", "Error:", "Donor weights:") %in% out))
  expect_false("Time variation of coefficients:" %in% out)
  labels <- sub("^[0-9]+ +([^ ]+).*", "\\1", out)
  expect_true(all(c("intercept", "x", "sigma") %in% labels))
  expect_true(any(grepl("Effective number of donors", out, fixed = TRUE)))
})

test_that("print.bscmfit shows the convergence measures unless disabled", {
  suppressWarnings(out <- capture.output(print(fit1_x)))
  expect_true(any(grepl("mcse_mean", out, fixed = TRUE)))
  suppressWarnings(
    out <- capture.output(print(fit1_x, diagnostics = FALSE))
  )
  expect_false(any(grepl("mcse_mean", out, fixed = TRUE)))
})

test_that("print.bscmfit shows the treatment effects and the model fit", {
  suppressWarnings(out <- capture.output(print(fit1_x)))
  expect_true(
    any(grepl("Average post-treatment effect", out, fixed = TRUE))
  )
  expect_true(any(grepl("Bayesian R2", out, fixed = TRUE)))
  expect_true(any(grepl("Post-treatment RMSE", out, fixed = TRUE)))
})

test_that("print.bscmfit omits the unit-specific statistics for N > 1", {
  suppressWarnings(out <- capture.output(print(fitN_xz)))
  expect_true("Coefficients:" %in% out)
  expect_false("Error:" %in% out)
  expect_true("Donor weights:" %in% out)
  expect_true(
    any(grepl("Average effective number of donors", out, fixed = TRUE))
  )
  expect_false(any(grepl("intercept", out, fixed = TRUE)))
  expect_true(
    any(grepl("available via `summary()`", out, fixed = TRUE))
  )
  expect_true(
    any(grepl("Average post-treatment effect", out, fixed = TRUE))
  )
  expect_false(any(grepl("Bayesian R2", out, fixed = TRUE)))
})

test_that("print.bscmfit omits the treatment effects without the data", {
  x <- fit1_x
  x$data <- NULL
  suppressWarnings(out <- capture.output(print(x)))
  expect_true("Error:" %in% out)
  expect_true(any(grepl("Some of the summary statistics", out, fixed = TRUE)))
  expect_false(
    any(grepl("Average post-treatment effect", out, fixed = TRUE))
  )
})
