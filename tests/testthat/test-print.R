test_that("print.bscmfit returns x invisibly", {
  out <- capture.output(result <- print(fit1_int))
  expect_identical(result, fit1_int)
})

test_that("print.bscmfit output contains expected header", {
  out <- capture.output(print(fit1_int))
  expect_true(any(grepl("Bayesian synthetic control model", out)))
})

test_that("print.bscmfit output contains MCMC info", {
  out <- capture.output(print(fit1_int))
  expect_true(any(grepl("MCMC sampling", out)))
})

test_that("print.bscmfit shows treated unit info for single treated", {
  out <- capture.output(print(fit1_int))
  expect_true(any(grepl("Treated unit:", out)))
  expect_true(any(grepl("Number of donors:", out)))
})

test_that("print.bscmfit shows number of treated units for N > 1", {
  suppressWarnings(out <- capture.output(print(fitN_int)))
  expect_true(any(grepl("Number of treated units:", out)))
})

test_that("print.bscmfit shows AR(1) residuals for AR(1) model", {
  suppressWarnings(out <- capture.output(print(fit1_noint)))
  expect_true(any(grepl("AR(1) residuals", out, fixed = TRUE)))
})

test_that("print.bscmfit does not show AR(1) for IID model", {
  out <- capture.output(print(fit1_int))
  expect_false(any(grepl("AR(1) residuals", out, fixed = TRUE)))
})

test_that("print.summary_bscmfit returns x invisibly", {
  suppressWarnings(s <- summary(fit1_int))
  out <- capture.output(result <- print(s))
  expect_identical(result, s)
})

test_that("print.summary_bscmfit produces output without error", {
  suppressWarnings(s <- summary(fit1_int))
  expect_output(print(s))
})
