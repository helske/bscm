test_that("coef() returns data.frame for intercept-only model", {
  d <- coef(fit_intercept)
  expect_s3_class(d, "draws_summary")
  expect_equal(d$variable, "Intercept")
})

test_that("coef() returns both intercept and beta for predictor model", {
  suppressWarnings(d <- coef(fit_x))
  expect_s3_class(d, "draws_summary")
  expect_equal(d$variable, c("Intercept", "Coef_x"))
})

test_that("coef() errors for no-intercept no-predictor model", {
  expect_error(
    coef(fit_no_intercept),
    "The model does not contain an intercept or any predictors"
  )
})

test_that("coef() validates probs argument", {
  expect_error(
    coef(fit_intercept, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
  expect_error(
    coef(fit_intercept, probs = 1.5),
    "Argument `probs` must be a <numeric> vector"
  )
  expect_error(
    coef(fit_intercept, probs = -0.1),
    "Argument `probs` must be a <numeric> vector with values between 0 and 1"
  )
})

test_that("coef() respects custom probs", {
  d <- coef(fit_intercept, probs = c(0.1, 0.9))
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
})
