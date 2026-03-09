test_that("as.data.frame returns a data.frame", {
  d <- as.data.frame(fit_intercept)
  expect_s3_class(d, "data.frame")
  expect_equal(dim(d), c(ndraws(fit_intercept), 32L))
  expect_equal(
    names(d), 
    c("alpha", "sigma", paste0("omega[", seq_len(30), "]"))
  )
  d <- as.data.frame(fit_x, parameters = c("beta", "lp__"))
  expect_equal(dim(d), c(ndraws(fit_intercept), 2L))
  expect_equal(names(d), c("beta[1]", "lp__"))
  expect_error(
    as.data.frame(fit_intercept, parameters = c("none", "aa", "something")),
    'Model does not contain parameters "none", "aa", and "something"\\.'
  )
})
