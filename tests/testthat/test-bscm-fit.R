test_that("bscm() returns a bscmfit object ", {
  expect_s3_class(fit1_int, "bscmfit")
  expected_names <- c(
    "stanfit",
    "data",
    "setup",
    "converge",
    "elapsed_time",
    "call"
  )
  expect_true(all(expected_names %in% names(fit1_int)))
})

test_that("bscmfit$setup has correct structure", {
  setup <- fit1_int$setup
  expect_type(setup, "list")
  expected_fields <- c(
    "outcome",
    "treatment",
    "treated",
    "donors",
    "unit",
    "time",
    "times",
    "T_pre",
    "T_total",
    "has_intercept",
    "predictors",
    "beta_names",
    "gamma_names",
    "model_type",
    "priors"
  )
  expect_true(all(expected_fields %in% names(setup)))
})

test_that("bscmfit$setup values are correct", {
  setup <- fit1_int$setup
  expect_equal(setup$outcome, "y")
  expect_equal(setup$treatment, "treatment")
  expect_equal(setup$treated, "1")
  expect_equal(setup$unit, "id")
  expect_equal(setup$time, "time")
  expect_equal(setup$T_pre, 30L)
  expect_equal(setup$T_total, 35L)
  expect_true(setup$has_intercept)
  expect_equal(setup$predictors, character(0))
  expect_null(setup$beta_names)
  expect_equal(length(setup$donors), 50L)
  expect_s3_class(fit1_int$data, "data.frame")
  expect_equal(get_predictors(fitN_xz), c("x", "z"))
  expect_true(has_intercept(fitN_xz))
  expect_equal(fit1_int$setup$model_type, "bscm_a1_x0_w0_dr_iid")
  expect_equal(fit1_noint$setup$model_type, "bscm_a0_x0_w0_dr_iid")
  expect_equal(fit1_x$setup$model_type, "bscm_a1_x1_w0_dr_iid")
  expect_equal(fitN_xz$setup$model_type, "bscm_a1_x1_w0_dr_iid")
  expect_equal(get_treated(fitN_xz), c("1", "2", "3"))
})

test_that("bscmfit$setup with tv() has correct structure", {
  setup <- fitN_tv$setup
  expect_equal(setup$model_type, "bscm_a1_x1_w1_dr_iid")
  expect_identical(setup$gamma_names, "x")
  expect_true("x" %in% setup$beta_names)
  expect_true("z" %in% setup$beta_names)
})

test_that("bscmfit$call is recorded", {
  expect_true(is.call(fit1_int$call))
})

test_that("bscmfit$converge is a bscmfit_diagnostics object", {
  expect_s3_class(fit1_int$converge, "bscmfit_diagnostics")
})
