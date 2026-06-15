test_that("plot_residuals returns ggplot for single treated unit", {
  p <- plot_residuals(fit1_int)
  expect_s3_class(p, "ggplot")
})

test_that("plot_residuals returns named list of ggplots for N > 1", {
  suppressWarnings(plots <- plot_residuals(fitN_int))
  expect_type(plots, "list")
  expect_named(plots, get_treated(fitN_int))
  expect_true(all(vapply(plots, inherits, logical(1L), "ggplot")))
})

test_that("plot_residuals type = 'autocorrelation' returns ggplot", {
  p <- plot_residuals(fit1_int, type = "autocorrelation")
  expect_s3_class(p, "ggplot")
})

test_that("plot_residuals validates type argument", {
  expect_error(
    plot_residuals(fit1_int, type = "invalid"),
    "Argument `type` must be either"
  )
})

test_that("plot_residuals validates probs argument length", {
  expect_error(
    plot_residuals(fit1_int, probs = c(0.025, 0.5, 0.975)),
    "Argument `probs` must be a"
  )
})

test_that("residual_acf returns data.frame with correct structure", {
  d <- residual_acf(fit1_int, max_lag = 4L)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 4L)
  expect_true("lag" %in% names(d))
  expect_true("mean" %in% names(d))
  expect_equal(d$lag, seq_len(4L))
})

test_that("residual_acf accepts draws directly and matches bscmfit result", {
  d_raw <- residuals(fit1_int, summary = FALSE)
  d_from_draws <- residual_acf(d_raw, max_lag = 3L)
  d_from_fit <- residual_acf(fit1_int, max_lag = 3L)
  expect_equal(d_from_draws, d_from_fit)
})

test_that("residual_acf validates max_lag argument", {
  expect_error(
    residual_acf(fit1_int, max_lag = -1L),
    "Argument `max_lag` must be a single positive integer"
  )
  expect_error(
    residual_acf(fit1_int, max_lag = 1.5),
    "Argument `max_lag` must be a single positive integer"
  )
})
