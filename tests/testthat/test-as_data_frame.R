test_that("as.data.frame returns a data.frame", {
  d <- as.data.frame(fit1_int)
  expect_s3_class(d, "data.frame")
  expect_equal(dim(d), c(ndraws(fit1_int), 32L))
  expect_equal(
    names(d), 
    c("alpha[1]", "sigma[1]", paste0("omega[1,", seq_len(30), "]"))
  )
  d <- as.data.frame(fit1_x, parameters = c("beta", "lp__"))
  expect_equal(dim(d), c(ndraws(fit1_int), 2L))
  expect_equal(names(d), c("beta[1]", "lp__"))
  
  d <- as.data.frame(
    fitN_xz, parameters = c("beta", "y_rep")
  )
  expect_equal(dim(d), c(ndraws(fitN_xz), 122L))
  expect_equal(
    names(d), 
    c(paste0("beta[", 1:2, "]"),
      paste0("y_rep[", seq_len(40), ",", rep(1:3, each = 40), "]")
    )
  )
  
  expect_error(
    as.data.frame(fit1_int, parameters = c("sigma", "none", "aa", "something")),
    'Model does not contain parameters "none", "aa", and "something"\\.'
  )
})
