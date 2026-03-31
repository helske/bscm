test_that("covariate_imbalance has expected columns", {
  d <- covariate_imbalance(fit1_x, plot = FALSE)
  expect_s3_class(d, "draws_summary")
  expect_named(d, c("id", "time", "mean", "q2.5", "q97.5"))
  d <- covariate_imbalance(fitN_xz, plot = FALSE)
  expect_s3_class(d, "draws_summary")
  expect_named(d, c("id", "time", "mean", "q2.5", "q97.5"))
  expect_equal(nrow(d), 120L)
})

test_that("covariate_imbalance errors for model without predictors", {
  expect_error(
    covariate_imbalance(fit1_int, plot = FALSE),
    "The model does not contain any predictors"
  )
})

test_that("covariate_imbalance validates plot argument", {
  expect_error(
    covariate_imbalance(fit1_x, plot = "yes"),
    "Argument `plot` must be a single <logical> value"
  )
})

test_that("covariate_imbalance validates probs argument", {
  expect_error(
    covariate_imbalance(fit1_x, plot = FALSE, probs = "x"),
    "Argument `probs` must be a <numeric> vector"
  )
})
