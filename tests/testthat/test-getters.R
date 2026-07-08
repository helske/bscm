# Tests for get_stanfit, ndraws, nchains

test_that("get_stanfit returns a stanfit object", {
  sf <- get_stanfit(fit1_int)
  expect_s4_class(sf, "stanfit")
})

test_that("get_stanfit.default extracts stanfit from list", {
  x <- list(stanfit = fit1_int$stanfit)
  sf <- get_stanfit(x)
  expect_s4_class(sf, "stanfit")
})

test_that("ndraws returns a positive integer", {
  nd <- ndraws(fit1_int)
  expect_identical(nd, 2000L)
})

test_that("nchains returns the correct integer", {
  nc <- nchains(fit1_int)
  expect_identical(nc, 2L)
})

test_that("get_treated returns character vector of treated unit IDs", {
  expect_type(get_treated(fit1_int), "character")
  expect_length(get_treated(fit1_int), 1L)
  expect_length(get_treated(fitN_int), 3L)
})

test_that("get_donors returns a character vector", {
  donors <- get_donors(fit1_int)
  expect_type(donors, "character")
  expect_length(donors, get_J(fit1_int))
})

test_that("get_times length matches get_T_total", {
  expect_equal(length(get_times(fit1_int)), get_T_total(fit1_int))
})

test_that("get_T_pre returns a named integer", {
  tp <- get_T_pre(fit1_int)
  expect_type(tp, "integer")
  expect_named(tp, get_treated(fit1_int))
})

test_that("get_N returns 1 for single treated and 3 for multiple treated", {
  expect_equal(get_N(fit1_int), 1L)
  expect_equal(get_N(fitN_int), 3L)
})

test_that("get_J returns a positive integer", {
  j <- get_J(fit1_int)
  expect_type(j, "integer")
  expect_gt(j, 0L)
})

test_that("get_outcome returns outcome variable name from bscmfit", {
  expect_equal(get_outcome(fit1_int), "y")
})

test_that("get_outcome works on a formula", {
  expect_equal(get_outcome(y ~ x + z), "y")
})

test_that("get_outcome errors on non-formula input", {
  expect_error(get_outcome("y ~ x"), "must be a <formula> object")
})

test_that("get_outcome errors on one-sided formula", {
  expect_error(get_outcome(~x), "with an outcome variable")
})

test_that("get_outcome errors on multi-variable LHS", {
  expect_error(get_outcome(cbind(a, b) ~ x), "one outcome variable")
})

test_that("has_intercept is TRUE for intercept models and FALSE otherwise", {
  expect_true(has_intercept(fit1_int))
  expect_false(has_intercept(fit1_noint))
})

test_that("has_predictors is FALSE for intercept-only and TRUE for covariate models", {
  expect_false(has_predictors(fit1_int))
  expect_true(has_predictors(fit1_x))
})

test_that("has_tv_coefs is FALSE for standard models and TRUE for tv() models", {
  expect_false(has_tv_coefs(fit1_int))
  expect_true(has_tv_coefs(fitN_tv))
})


test_that("get_standata returns a named list with expected Stan data fields", {
  x <- get_standata(fit1_int)
  expect_type(x, "list")
  expect_true(all(c("J", "N", "T", "T_pre", "y", "Z") %in% names(x)))
  expect_equal(x$N, get_N(fit1_int))
  expect_equal(x$J, get_J(fit1_int))
})

test_that("get_standata errors when no saved data", {
  fit_nosave <- fit1_int
  fit_nosave$data <- NULL
  expect_error(get_standata(fit_nosave), "save_data = TRUE")
})

