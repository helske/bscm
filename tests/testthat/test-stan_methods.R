test_that("posterior_predict returns a matrix", {
  pp <- posterior_predict(fit1_int)
  expect_true(is.matrix(pp))
})

test_that("posterior_predict has correct dimensions", {
  pp <- posterior_predict(fit1_int)
  expect_equal(dim(pp), c(ndraws(fit1_int), get_T_total(fit1_int)))
  pp <- posterior_predict(fitN_int)
  expect_equal(dim(pp), c(ndraws(fitN_int), 3 * get_T_total(fitN_int)))
})

test_that("posterior_epred returns a matrix with correct dimensions", {
  pe <- posterior_epred(fit1_x)
  expect_true(is.matrix(pe))
  expect_equal(dim(pe), c(ndraws(fit1_x), get_T_total(fit1_x)))
  expect_equal(
    colnames(pe),
    paste0("y_mean[", seq_len(get_T_total(fit1_x)), ",1]")
  )
})

test_that("posterior_linpred returns same as posterior_epred", {
  pe <- posterior_epred(fitN_xz)
  pl <- posterior_linpred(fitN_xz)
  expect_identical(pl, pe)
})

test_that("posterior_linpred ignores transform argument", {
  pl1 <- posterior_linpred(fit1_int, transform = FALSE)
  pl2 <- posterior_linpred(fit1_int, transform = TRUE)
  expect_equal(pl1, pl2)
})

test_that("log_lik returns a matrix", {
  ll <- log_lik(fit1_int)
  expect_true(is.matrix(ll))
})

test_that("log_lik has correct dimensions", {
  ll <- log_lik(fit1_int)
  expect_equal(dim(ll), unname(c(ndraws(fit1_int), get_T_pre(fit1_int))))
  ll <- log_lik(fitN_int)
  expect_equal(dim(ll), c(ndraws(fitN_int), sum(get_T_pre(fitN_int))))
})
