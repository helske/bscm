test_that("covariate_distance has expected columns", {
  suppressWarnings(d <- covariate_distance(fit_x, plot = FALSE))
  expect_s3_class(d, "draws_summary")
})

test_that("covariate_distance errors for model without predictors", {
  expect_error(
    covariate_distance(fit_intercept, plot = FALSE),
    "The model does not contain any predictors"
  )
})

test_that("covariate_distance validates plot argument", {
  expect_error(
    covariate_distance(fit_x, plot = "yes"),
    "Argument `plot` must be a single <logical> value"
  )
})

test_that("covariate_distance validates probs argument", {
  expect_error(
    covariate_distance(fit_x, plot = FALSE, probs = "x"),
    "Argument `probs` must be a <numeric> vector"
  )
})
