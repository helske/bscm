test_that("fitted returns data.frame with correct structure", {
  suppressWarnings(d <- fitted(fit1_noint))
  expect_s3_class(d, "tbl_df")
  expect_identical(names(d)[1:4], c("time", "treatment", "mean", "sd"))
  expect_equal(nrow(d), get_T_total(fit1_noint))

  suppressWarnings(d <- fitted(fitN_int))
  expect_equal(d$time, rep(get_times(fitN_int), 3))
  expect_equal(d$id, rep(get_treated(fitN_int), each = get_T_total(fitN_int)))
})

test_that("fitted returns posterior samples when summary = FALSE", {
  d <- fitted(fit1_int, summary = FALSE)
  expect_s3_class(d, "tbl_df")
  expect_true(posterior::is_rvar(d$y_mean))
  expect_equal(nrow(d), get_T_total(fit1_int))
})

test_that("fitted respects custom probs", {
  d <- fitted(fit1_int, probs = c(0.1, 0.9))
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
})

test_that("fitted validates summary argument", {
  expect_error(
    fitted(fit1_int, summary = "yes"),
    "Argument `summary` must be a single"
  )
})

test_that("fitted validates probs argument", {
  expect_error(
    fitted(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("residuals equal observed minus fitted", {
  f <- fitted(fit1_int, summary = FALSE)
  r <- residuals(fit1_int, summary = FALSE, pretreatment_only = FALSE)
  y <- c(get_stan_y(fit1_int))
  expect_equal(
    mean(posterior::E(r$residuals - (y - f$y_mean))),
    0,
    tolerance = 1e-10
  )
})
