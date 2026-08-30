test_that("as_draws returns a draws_array", {
  d <- as_draws(fit1_int, "alpha")
  expect_s3_class(d, "draws_array")
  expect_equal(dim(d), c(ndraws(fit1_int) / 2L, 2L, 1L))
  expect_equal(posterior::variables(d), "alpha[1]")
  d <- as_draws(fit1_int, "omega", include = FALSE)
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

test_that("as_draws_rvars returns a draws_rvars", {
  d <- as_draws_rvars(fit1_int, "alpha")
  expect_s3_class(d, "draws_rvars")
  expect_equal(posterior::variables(d), "alpha")
  expect_equal(dim(d$alpha), 1L)
  expect_equal(ndraws(d), ndraws(fit1_int))
  expect_equal(nchains(d), 2L)
  d <- as_draws_rvars(fit1_int, "omega", include = FALSE)
  expect_equal(posterior::variables(d), c("sigma", "alpha", "lp__"))
})
