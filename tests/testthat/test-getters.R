# Tests for get_stanfit, ndraws, nchains

test_that("get_stanfit returns a stanfit object", {
  sf <- get_stanfit(fit_intercept)
  expect_s4_class(sf, "stanfit")
})

test_that("get_stanfit.default extracts stanfit from list", {
  x <- list(stanfit = fit_intercept$stanfit)
  sf <- get_stanfit(x)
  expect_s4_class(sf, "stanfit")
})

test_that("ndraws returns a positive integer", {
  nd <- ndraws(fit_intercept)
  expect_type(nd, "integer")
  expect_equal(nd, 150L)
})

test_that("nchains returns the correct integer", {
  nc <- nchains(fit_intercept)
  expect_type(nc, "integer")
  expect_equal(nc, 1L)
})
