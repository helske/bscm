test_that("plot.bscmfit returns a ggplot object", {
  p <- plot(fit1_int)
  expect_s3_class(p, "ggplot")
})

test_that("plot.bscmfit returns averaged ggplot for N > 1", {
  p <- plot(fitN_int)
  expect_s3_class(p[["2"]], "ggplot")
})

test_that("plot.bscmfit works with custom probs", {
  p <- plot(fit1_int, probs = c(0.1, 0.9))
  expect_s3_class(p, "ggplot")
})

test_that("plot.bscmfit validates probs argument", {
  expect_error(
    plot(fit1_int, probs = c(0.025)),
    "Argument `probs` must be a <numeric> vector of length 2"
  )
  expect_error(
    plot(fit1_int, probs = c(-0.1, 0.975)),
    "Argument `probs` must be a <numeric> vector of length 2"
  )
  expect_error(
    plot(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector of length 2"
  )
})

test_that("plot_coefs returns a ggplot for both coefficient types", {
  expect_s3_class(plot_coefs(fit1_x), "ggplot")
  suppressWarnings(expect_s3_class(plot_coefs(fitN_tv), "ggplot"))
  suppressWarnings(
    expect_s3_class(plot_coefs(fitN_tv, type = "fixed"), "ggplot")
  )
})

test_that("plot_coefs errors for a model without predictors", {
  expect_error(
    plot_coefs(fit1_int),
    "The model does not contain any predictors"
  )
})

test_that("plot_coefs does not show the intercepts", {
  d <- ggplot2::ggplot_build(plot_coefs(fit1_x))$plot$data
  expect_false("Intercept" %in% as.character(d$parameter))
})

test_that("plot_coefs errors for time-varying coefficients of a model with none", {
  expect_error(
    plot_coefs(fit1_x, type = "varying"),
    "The model does not contain any time-varying predictors"
  )
})
