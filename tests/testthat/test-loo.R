test_that("loo returns a loo object", {
  l <- suppressWarnings(loo(fit1_int))
  expect_s3_class(l, "loo")
  expect_true("estimates" %in% names(l))
  expect_true("pointwise" %in% names(l))
  d <- attr(l, "diagnostics")
  expect_s3_class(d, "tbl_df")
  expect_named(d, c("id", "time", "pareto_k", "n_eff", "r_eff"))
})

test_that("loo works with r_eff = FALSE", {
  l <- suppressWarnings(loo(fit1_int, r_eff = FALSE))
  expect_s3_class(l, "loo")
})

test_that("loo works for multiple treated units", {
  l <- suppressWarnings(loo(fitN_tv))
  expect_s3_class(l, "loo")
  expect_true("estimates" %in% names(l))
})
