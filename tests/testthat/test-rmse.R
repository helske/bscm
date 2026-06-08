test_that("rmse returns data.frame with RMSEs", {
  d <- rmse(fit1_int)
  expect_s3_class(d, "tbl_df")
  expect_equal(d$treatment, 0:1)

  d <- rmse(fitN_xz, probs = numeric(0))
  expect_s3_class(d, "tbl_df")
  expect_named(
    d,
    c(
      "id",
      "treatment",
      "mean",
      "sd",
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean"
    )
  )
  d <- rmse(fitN_xz, probs = 0.5, average = TRUE)
  expect_s3_class(d, "tbl_df")
  expect_named(
    d,
    c(
      "treatment",
      "mean",
      "sd",
      "q50",
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean"
    )
  )
})

test_that("rmse respects custom probs", {
  d <- rmse(fit1_int, probs = c(0.1, 0.9))
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
})

test_that("rmse validates probs argument", {
  expect_error(
    rmse(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("rmse average argument works", {
  expect_error(
    rmse(fit1_int, average = 4),
    "Argument `average` must be a single"
  )
  d1 <- rmse(fit1_int)
  d2 <- rmse(fit1_int, average = TRUE)
  expect_equal(d1, d2)
  d1 <- rmse(fitN_int)
  d2 <- rmse(fitN_int, average = FALSE)
  expect_equal(d1, d2)
  d2 <- rmse(fitN_int, average = TRUE)
  expect_identical(dim(d2), c(2L, 9L))
})
