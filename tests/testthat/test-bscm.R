test_that("Incorrect arguments to bscm result in meaningful error message", {
  expect_error(
    bscm(),
    "Argument `formula` is missing\\."
  )
  expect_error(
    bscm(r ~ 1, data = simulated_data, treatment = "treatment"),
    "Can't find outcome variable `r` in `data`\\."
  )
  expect_error(
    bscm( ~ 1, data = simulated_data, treatment = "treatment"),
    "Argument `formula` must be a <formula> object with an outcome variable on the left-hand side\\."
  )
  expect_error(
    bscm(c(y, x) ~ 1, data = simulated_data, treatment = "treatment"),
    "Argument `formula` must be a <formula> object with one outcome variable on the left-hand side\\."
  )
  expect_error(
    bscm(y ~ 1),
    "Argument `data` is missing\\."
  )
  expect_error(
    bscm(y ~ 1, data = 5),
    "Argument `data` must be a <data\\.frame> object\\."
  )
  expect_error(
    bscm(y ~ 1, data = simulated_data),
    "Argument `treatment` is missing\\."
  )
  expect_error(
    bscm(y ~ 1, data = simulated_data, treatment = "z"),
    "Can't find treatment variable `z` in `data`\\."
  )
  expect_error(
    bscm(y ~ 1, data = simulated_data, treatment = "x"),
    "Variable `x` in `data` should contain either logical or binary values indicating the pre- and post-treatment time points\\."
  )
  expect_error(
    bscm(y ~ 1, data = simulated_data[, -1], treatment = "treatment"),
    "Can't find time index variable `time` in `data`\\."
  )
  expect_error(
    bscm(y ~ 1, data = simulated_data, treatment = "treatment", time = "t"),
    "Can't find time index variable `t` in `data`\\."
  )
  expect_error(
    bscm(y ~ 1, data = simulated_data[, -2], treatment = "treatment"),
    "Can't find unit index variable `id` in `data`\\."
  )
})
