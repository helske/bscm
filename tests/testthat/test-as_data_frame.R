test_that("as.data.frame returns a data.frame", {
  d <- as.data.frame(fit1_int)
  J <- get_J(fit1_int)
  expect_s3_class(d, "data.frame")
  expect_equal(dim(d), c(ndraws(fit1_int), 2))
  expect_equal(names(d), c("alpha[1]", "sigma[1]"))
  d <- as.data.frame(fit1_x, parameters = c("beta", "lp__"))
  expect_equal(dim(d), c(ndraws(fit1_x), 2L))
  expect_equal(names(d), c("beta[1]", "lp__"))

  d <- as.data.frame(
    fitN_xz,
    parameters = c("beta", "sigma[1]")
  )
  expect_equal(dim(d), c(ndraws(fitN_xz), 3L))
  expect_equal(names(d), c("beta[1]", "beta[2]", "sigma[1]"))

  expect_error(
    as.data.frame(fit1_int, parameters = c("sigma", "none", "aa", "something")),
    'Model does not contain parameters "none", "aa", and "something"\\.'
  )
})
