test_that("summary returns a bscmfit_summary", {
  suppressWarnings(s <- summary(fit1_int))
  expect_s3_class(s, "summary_bscmfit")
  expect_equal(nrow(s), 6L)
})

test_that("summary validates probs argument", {
  expect_error(
    summary(fit1_noint, probs = "x"),
    "Argument `probs` must be a <numeric> vector"
  )
})
