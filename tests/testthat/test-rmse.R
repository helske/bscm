test_that("rmse returns data.frame with RMSEs", {
  d <- rmse(fit1_int)
  expect_s3_class(d, "data.frame")
  expect_equal(d$variable, c("Pre-RMSE", "Post-RMSE", "RMSE ratio"))

  d <- rmse(fitN_xz, probs = numeric(0))
  expect_s3_class(d, "data.frame")
  expect_named(
    d,
    c("variable", "mean", "sd", "rhat", "ess_bulk", "ess_tail")
  )
  d <- rmse(fitN_xz, probs = 0.5, average = FALSE)
  expect_s3_class(d, "data.frame")
  expect_named(
    d,
    c("id", "variable", "mean", "sd", "q50", "rhat", "ess_bulk", "ess_tail")
  )
  expect_equal(d$id, as.character(rep(1:3, times = 3)))
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
  d2 <- rmse(fitN_int, average = TRUE)
  expect_equal(d1, d2)
  d2 <- rmse(fitN_int, average = FALSE)
  expect_identical(dim(d2), c(9L, 9L))
})
