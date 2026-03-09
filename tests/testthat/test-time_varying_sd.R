test_that("time_varying_sd errors for non-time-varying model", {
  expect_error(
    time_varying_sd(fit_x, plot = FALSE),
    "The model was not estimated with `time_varying_effects = TRUE`"
  )
  expect_error(
    time_varying_sd(fit_intercept, plot = FALSE),
    "The model was not estimated with `time_varying_effects = TRUE`"
  )
})

test_that("time_varying_sd validates plot argument", {
  expect_error(
    time_varying_sd(fit_tv_x, plot = "yes"),
    "Argument `plot` must be a single <logical> value"
  )
})

test_that("time_varying_sd validates probs argument", {
  expect_error(
    time_varying_sd(fit_tv_x, plot = FALSE, probs = "x"),
    "Argument `probs` must be a <numeric> vector"
  )
})
test_that("time_varying_sd works", {
  expect_error(
    out <- time_varying_sd(fit_tv_x, plot = FALSE, probs = c(0.025, 0.9)),
    NA
  )
  expect_s3_class(out, "data.frame")
  expect_equal(dim(out), c(50L, 8L))
})
