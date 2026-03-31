test_that("sigma returns data.frame with correct structure", {
  d <- sigma(fit1_int)
  expect_s3_class(d, "data.frame")
  expect_equal(d$variable, "sigma")
  expect_true("mean" %in% names(d))
  expect_true("sd" %in% names(d))
  expect_equal(nrow(d), 1L)
})

test_that("sigma respects custom probs", {
  d <- sigma(fit1_int, probs = c(0.1, 0.9))
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
})

test_that("sigma validates probs argument", {
  expect_error(
    sigma(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("sigma returns positive values", {
  d <- sigma(fit1_int)
  expect_gt(d$q2.5, 0)
})
