test_that("residuals returns data.frame with correct structure", {
  suppressWarnings(d <- residuals(fit1_noint))
  expect_s3_class(d, "tbl_df")
  expect_identical(names(d)[1:4], c("time", "treatment", "mean", "sd"))
  expect_equal(nrow(d), unname(get_T_pre(fit1_noint)[1]))

  suppressWarnings(d <- residuals(fitN_int))
  expect_equal(
    d$id,
    rep(get_treated(fitN_int), each = get_T_pre(fitN_int)[1])
  )
})

test_that("residuals returns posterior draws when summary = FALSE", {
  d <- residuals(fit1_int, summary = FALSE)
  expect_s3_class(d, "tbl_df")
  expect_true(posterior::is_rvar(d$residuals))
  expect_equal(nrow(d), unname(get_T_pre(fit1_int)[1]))
})

test_that("residuals respects custom probs", {
  d <- residuals(fit1_int, probs = c(0.1, 0.9))
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
})

test_that("residuals validates summary argument", {
  expect_error(
    residuals(fit1_int, summary = "yes"),
    "Argument `summary` must be a single"
  )
})

test_that("residuals validates probs argument", {
  expect_error(
    residuals(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("residuals are reasonable", {
  suppressWarnings(d <- residuals(fit1_x))
  pre_resid <- d$mean[d$treatment == 0]
  expect_lt(max(abs(pre_resid)), 1)
})
