test_that("plot_weights returns ggplot for single treated unit", {
  p <- plot_weights(fit1_int)
  expect_s3_class(p, "ggplot")
})

test_that("plot_weights returns named list of ggplots for N > 1", {
  suppressWarnings(plots <- plot_weights(fitN_int))
  expect_type(plots, "list")
  expect_named(plots, get_treated(fitN_int))
  expect_true(all(vapply(plots, inherits, logical(1L), "ggplot")))
})

test_that("plot_weights works with point_estimate = 'mean'", {
  p <- plot_weights(fit1_int, point_estimate = "mean")
  expect_s3_class(p, "ggplot")
})

test_that("plot_weights works with order options", {
  p <- plot_weights(fit1_int, order = "ascending")
  expect_s3_class(p, "ggplot")
  p <- plot_weights(fit1_int, order = "descending")
  expect_s3_class(p, "ggplot")
})

test_that("plot_weights validates point_estimate argument", {
  expect_error(
    plot_weights(fit1_int, point_estimate = "mode"),
    "Argument `point_estimate` must be either"
  )
})

test_that("plot_weights validates coverage argument", {
  expect_error(
    plot_weights(fit1_int, coverage = c(0, 0.95)),
    "Argument `coverage` must have values between 0 and 1"
  )
  expect_error(
    plot_weights(fit1_int, coverage = 1.1),
    "Argument `coverage` must have values between 0 and 1"
  )
})

test_that("plot_weights validates linewidth argument", {
  expect_error(
    plot_weights(fit1_int, linewidth = -1),
    "Argument `linewidth` must be a single positive number"
  )
})

test_that("plot_weights validates point_size argument", {
  expect_error(
    plot_weights(fit1_int, point_size = -1),
    "Argument `point_size` must be a single positive number"
  )
})

test_that("plot_weights validates reverse argument", {
  expect_error(
    plot_weights(fit1_int, reverse = "yes"),
    "Argument `reverse`"
  )
})
