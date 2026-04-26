test_that("plot_effects.bscmfit returns a ggplot for single treated unit", {
  p <- plot_effects(fit1_int)
  expect_s3_class(p, "ggplot")
})

test_that("plot_effects.bscmfit returns a list of ggplots for multiple treated units", {
  plots <- plot_effects(fitN_int)
  expect_type(plots, "list")
  expect_true(all(vapply(plots, inherits, logical(1L), "ggplot")))
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

test_that("plot_effects.bscm_ldo returns a ggplot with LDO lines", {
  ldo <- suppressWarnings(
    leave_donor_out(fit1_int, refresh = 0, chains = 1, iter = 200, seed = 1)
  )
  p <- plot_effects(ldo)
  expect_s3_class(p, "ggplot")
  ldo_layer_data <- vapply(
    p$layers, \(l) inherits(l$geom, "GeomLine"), logical(1L)
  )
  expect_true(sum(ldo_layer_data) >= 2L)
})

test_that("plot_effects.bscm_ldo validates probs argument", {
  ldo <- suppressWarnings(
    leave_donor_out(fit1_int, refresh = 0, chains = 2, iter = 200, seed = 1)
  )
  expect_error(
    plot_effects(ldo, probs = c(0.5, 1.5)),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("plot_effects.bscm_ldo errors when probs not in stored quantiles", {
  ldo <- suppressWarnings(
    leave_donor_out(
      fit1_int, probs = c(0.1, 0.9),
      refresh = 0, chains = 2, iter = 200, seed = 1
    )
  )
  expect_error(
    plot_effects(ldo, probs = c(0.025, 0.975)),
    "Quantile columns matching"
  )
})
