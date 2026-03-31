test_that("loo returns a loo object", {
  l <- suppressWarnings(loo(fit1_int))
  expect_s3_class(l, "loo")
  expect_true("estimates" %in% names(l))
  expect_true("pointwise" %in% names(l))
})
