test_that("coef() returns data.frame for intercept-only model", {
  d <- coef(fit1_int)
  expect_s3_class(d, "tbl_df")
  expect_equal(d$variable, "Intercept")
})

test_that("coef() returns both intercept and beta for model with covariates", {
  suppressWarnings(d <- coef(fitN_xz, type = "beta"))
  expect_s3_class(d, "tbl_df")
  expect_equal(d$variable, c("beta_x", "beta_z"))
  suppressWarnings(d <- coef(fitN_xz))
  expect_s3_class(d, "tbl_df")
  expect_equal(
    d$variable,
    c("Intercept", "Intercept", "Intercept", "beta_x", "beta_z")
  )
})

test_that("coef() returns rho for AR(1) no-intercept model", {
  d <- coef(fit1_noint)
  expect_s3_class(d, "tbl_df")
  expect_equal(d$variable, "rho")
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
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
})
