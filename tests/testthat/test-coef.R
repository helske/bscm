test_that("coef() returns data.frame for intercept-only model", {
  d <- coef(fit1_int)
  expect_s3_class(d$alpha, "draws_summary")
  expect_equal(d$alpha$variable, "Intercept")
})

test_that("coef() returns both intercept and beta for predictor model", {
  suppressWarnings(d <- coef(fitN_xz))
  expect_s3_class(d$alpha, "draws_summary")
  expect_s3_class(d$beta, "draws_summary")
  expect_equal(d$beta$variable, c("Coef_x", "Coef_z"))
})

test_that("coef() errors for no-intercept no-predictor model", {
  expect_error(
    coef(fit1_noint),
    "The model does not contain an intercept or any predictors"
  )
})

test_that("coef() validates probs argument", {
  expect_error(
    coef(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
  expect_error(
    coef(fit1_int, probs = 1.5),
    "Argument `probs` must be a <numeric> vector"
  )
  expect_error(
    coef(fit1_int, probs = -0.1),
    "Argument `probs` must be a <numeric> vector with values between 0 and 1"
  )
})

test_that("coef() respects custom probs", {
  d <- coef(fit1_int, probs = c(0.1, 0.9))
  expect_true("q10" %in% names(d$alpha))
  expect_true("q90" %in% names(d$alpha))
})
