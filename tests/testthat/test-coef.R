test_that("coef() returns the intercept for an intercept-only model", {
  d <- coef(fit1_int)
  expect_s3_class(d, "tbl_df")
  expect_equal(d$variable, "Intercept")
})

test_that("coef() returns the constant coefficients by default", {
  suppressWarnings(d <- coef(fitN_xz))
  expect_s3_class(d, "tbl_df")
  expect_equal(d$variable, c("coef_x", "coef_z"))
})

test_that("coef() returns the intercepts of every treated unit", {
  suppressWarnings(d <- coef(fitN_xz, type = "intercept"))
  expect_equal(d$variable, rep("Intercept", get_N(fitN_xz)))
  expect_equal(d[[get_unit(fitN_xz)]], get_treated(fitN_xz))
})

test_that("coef() defaults to the fixed coefficients", {
  suppressWarnings(d <- coef(fitN_tv))
  expect_equal(d$variable, "coef_z")
})

test_that("coef() returns the total effect of a time-varying predictor", {
  suppressWarnings(d <- coef(fitN_tv, type = "varying"))
  expect_equal(unique(d$variable), "coef_x")
  expect_equal(nrow(d), get_T_total(fitN_tv))
  expect_true(get_time(fitN_tv) %in% names(d))
})

test_that("coef(type = 'varying') equals beta + gamma", {
  d <- coef(fitN_tv, type = "varying", summary = FALSE)
  beta <- as_draws_rvars(fitN_tv, "beta")$beta[fitN_tv$setup$tv_idx]
  gamma <- as_draws_rvars(fitN_tv, "gamma")$gamma
  total <- c(gamma) + rep(beta, each = get_T_total(fitN_tv))
  expect_equal(
    as.vector(posterior::draws_of(d$coef)),
    as.vector(posterior::draws_of(total))
  )
})

test_that("coef() errors for a model without an intercept and predictors", {
  expect_error(
    coef(fit1_noint),
    "The model does not contain an intercept nor any predictors"
  )
})

test_that("coef() errors for time-varying coefficients of a model with none", {
  expect_error(
    coef(fit1_x, type = "varying"),
    "The model does not contain any time-varying predictors"
  )
})

test_that("coef() validates type argument", {
  expect_error(
    coef(fit1_int, type = "beta"),
    'Argument `type` must be one of "intercept", "fixed", or "varying"'
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
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
})

test_that("coef() returns the draws as a single tibble", {
  d <- coef(fit1_x, summary = FALSE)
  expect_s3_class(d, "tbl_df")
  expect_named(d, c("parameter", "coef"))
  expect_equal(d$parameter, "coef_x")
  expect_true(posterior::is_rvar(d$coef))
  d <- coef(fit1_x, type = "intercept", summary = FALSE)
  expect_named(d, c("parameter", "id", "coef"))
})
