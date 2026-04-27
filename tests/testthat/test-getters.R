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
