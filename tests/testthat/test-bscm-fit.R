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
    "has_icpt",
    "has_x",
    "has_w",
    "has_ar1",
    "predictors",
    "beta_names",
    "gamma_names",
    "tv_idx",
    "spline_df",
    "knot_spacing",
    "prior_only",
    "marginalize_sigma",
    "excluded_pars"
  )
  expect_true(all(expected_fields %in% names(setup)))
})

test_that("sigma is marginalized by default", {
  expect_true(fit1_int$setup$marginalize_sigma)
  expect_true(fitN_tv$setup$marginalize_sigma)
  expect_equal(fit1_int$priors$sigma$distribution, "half_normal")
  sigma_draws <- as.matrix(get_stanfit(fitN_int), pars = "sigma")
  expect_equal(ncol(sigma_draws), 3L)
  expect_true(all(sigma_draws > 0))
})

test_that("incompatible prior for sigma turns off marginalization by default", {
  fit <- suppressWarnings(
    bscm(
      y ~ 1,
      data = single_treated,
      priors = list(sd = exponential_pr(1)),
      chains = 1,
      iter = 200,
      refresh = 0,
      seed = 1,
      mcmc_diagnostics = FALSE
    )
  )
  expect_false(fit$setup$marginalize_sigma)
  expect_error(
    bscm(
      y ~ 1,
      data = single_treated,
      priors = list(sigma = exponential_pr(1)),
      marginalize_sigma = TRUE,
      chains = 1,
      iter = 200,
      refresh = 0,
      seed = 1,
      mcmc_diagnostics = FALSE
    ),
    "is not supported with"
  )
})

test_that("bscmfit$setup values are correct", {
  setup <- fit1_int$setup
  expect_equal(setup$outcome, "y")
  expect_equal(setup$treatment, "treatment")
  expect_equal(setup$treated, "1")
  expect_equal(setup$unit, "id")
  expect_equal(setup$time, "time")
  expect_equal(setup$T_pre, 40L)
  expect_equal(setup$T_total, 50L)
  expect_true(setup$has_icpt)
  expect_equal(setup$predictors, character(0))
  expect_null(setup$beta_names)
  expect_equal(length(setup$donors), 50L)
  expect_s3_class(fit1_int$data, "data.frame")
  expect_equal(get_predictors(fitN_xz), c("x", "z"))
  expect_true(has_intercept(fitN_xz))
  expect_equal(get_treated(fitN_xz), c("1", "2", "3"))
  expect_equal(fit1_noint$setup$has_ar1, has_ar1_error(fit1_noint))
  expect_false(fit1_int$setup$has_ar1)
})

test_that("bscmfit$setup with tv() has correct structure", {
  setup <- fitN_tv$setup
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
