test_that("plot_effects.bscmfit returns a ggplot for single treated unit", {
  p <- plot_effects(fit1_int)
  expect_s3_class(p, "ggplot")
})

test_that("plot_effects.bscmfit returns averaged ggplot for N > 1", {
  p <- plot_effects(fitN_int)
  expect_s3_class(p, "ggplot")
})

test_that("plot_effects.bscmfit returns per-unit ggplot when unit is specified", {
  p <- plot_effects(fitN_int, unit = 1)
  expect_s3_class(p, "ggplot")
})

test_that("plot_effects.bscmfit validates unit argument", {
  expect_error(
    plot_effects(fitN_int, unit = 99),
    "Argument `unit` must be one of"
  )
})

test_that("plot_effects.bscmfit works with custom probs", {
  p <- plot_effects(fit1_int, probs = c(0.1, 0.9))
  expect_s3_class(p, "ggplot")
})

test_that("plot_effects.bscmfit validates probs argument", {
  expect_error(
    plot_effects(fit1_int, probs = 0.5),
    "Argument `probs` must be a <numeric> vector"
  )
  expect_error(
    plot_effects(fit1_int, probs = c(-0.1, 0.9)),
    "Argument `probs` must be a <numeric> vector"
  )
  expect_error(
    plot_effects(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
})
