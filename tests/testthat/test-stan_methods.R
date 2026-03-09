test_that("posterior_predict returns a matrix", {
  pp <- posterior_predict(fit_intercept)
  expect_true(is.matrix(pp))
})

test_that("posterior_predict has correct dimensions", {
  pp <- posterior_predict(fit_intercept)
  expect_equal(nrow(pp), ndraws(fit_intercept))
  expect_equal(ncol(pp), fit_intercept$setup$T_total)
})

test_that("posterior_epred returns a matrix with correct dimensions", {
  pe <- posterior_epred(fit_intercept)
  expect_true(is.matrix(pe))
  expect_equal(nrow(pe), ndraws(fit_intercept))
  expect_equal(ncol(pe), fit_intercept$setup$T_total)
})

test_that("posterior_linpred returns same as posterior_epred", {
  pe <- posterior_epred(fit_intercept)
  pl <- posterior_linpred(fit_intercept)
  expect_equal(pl, pe)
})

test_that("posterior_linpred ignores transform argument", {
  pl1 <- posterior_linpred(fit_intercept, transform = FALSE)
  pl2 <- posterior_linpred(fit_intercept, transform = TRUE)
  expect_equal(pl1, pl2)
})

test_that("log_lik returns a matrix", {
  ll <- log_lik(fit_intercept)
  expect_true(is.matrix(ll))
})

test_that("log_lik has correct dimensions (draws x T_pre)", {
  ll <- log_lik(fit_intercept)
  expect_equal(nrow(ll), ndraws(fit_intercept))
  expect_equal(ncol(ll), fit_intercept$setup$T_pre)
})
