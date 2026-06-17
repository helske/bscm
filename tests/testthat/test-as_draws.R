test_that("as_draws returns a draws_array", {
  d <- as_draws(fit1_int, "alpha")
  expect_s3_class(d, "draws_array")
  expect_equal(dim(d), c(ndraws(fit1_int) / 2L, 2L, 1L))
  expect_equal(posterior::variables(d),"alpha[1]")
  d <- as_draws(fit1_int, "omega",include = FALSE)
  expect_equal(
    posterior::variables(d),
    c("sigma[1]", "alpha[1]", "lp__")
  )
})

test_that("as_draws extracts rho for AR(1) model", {
  d <- as_draws(fit1_noint, "rho")
  expect_s3_class(d, "draws_array")
  expect_equal(dim(d), c(ndraws(fit1_noint) / 2L, 2L, 1L))
  expect_equal(posterior::variables(d), "rho[1]")
})
