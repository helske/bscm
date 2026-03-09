test_that("bscm() returns a bscmfit object ", {
  expect_s3_class(fit_intercept, "bscmfit")
  expected_names <- c(
    "stanfit", "data", "setup", "converge", "elapsed_time", "call"
  )
  expect_true(all(expected_names %in% names(fit_intercept)))
})

test_that("bscmfit$setup has correct structure", {
  setup <- fit_intercept$setup
  expect_type(setup, "list")
  expected_fields <- c(
    "outcome", "treatment", "treated", "donors", "unit", "time", "times",
    "T_pre", "T_total", "has_intercept", "predictors", "coef_names", "kappa",
    "model_type", "time_varying_effects", "priors"
  )
  expect_true(all(expected_fields %in% names(setup)))
})

test_that("bscmfit$setup values are correct for intercept-only model", {
  setup <- fit_intercept$setup
  expect_equal(setup$outcome, "y")
  expect_equal(setup$treatment, "treatment")
  expect_equal(setup$treated, "1")
  expect_equal(setup$unit, "id")
  expect_equal(setup$time, "time")
  expect_equal(setup$T_pre, 40L)
  expect_equal(setup$T_total, 50L)
  expect_true(setup$has_intercept)
  expect_equal(setup$predictors, character(0))
  expect_null(setup$coef_names)
  expect_false(setup$time_varying_effects)
  expect_equal(setup$kappa, 1)
  expect_equal(length(setup$donors), 30L)
  expect_s3_class(fit_intercept$data, "data.frame")
  expect_equal(setup$model_type, "bscm_int_nox_none")
})

test_that("bscmfit$setup values are correct for model with covariate", {
  setup <- fit_x$setup
  expect_equal(setup$predictors, "x")
  expect_equal(setup$coef_names, "x")
  expect_true(setup$has_intercept)
  expect_equal(setup$model_type, "bscm_int_x_const")
})

test_that("bscmfit$setup values for time-varying effect model", {
  setup <- fit_tv_x$setup
  expect_equal(setup$time_varying_effects, TRUE)
  expect_equal(setup$priors$pr_rate_tau, 1)
  expect_equal(setup$model_type, "bscm_int_x_varying")
})

test_that("bscmfit$call is recorded", {
  expect_true(is.call(fit_intercept$call))
})

test_that("bscmfit$converge is a bscmfit_diagnostics object", {
  expect_s3_class(fit_intercept$converge, "bscmfit_diagnostics")
})
