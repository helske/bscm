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
