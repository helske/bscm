test_that("placebo_effects validates arguments", {
  expect_error(
    placebo_effects(fit1_int, type = "invalid"),
    "Argument `type` must be either"
  )
  expect_error(
    placebo_effects(fit1_int, probs = "a"),
    "Argument `probs` must be a"
  )
  expect_error(
    placebo_effects(fit1_int, type = "time", L = 0),
    "Argument `L` must be a single integer between 2 and"
  )
  expect_error(
    placebo_effects(fit1_int, type = "time", L = "a"),
    "Argument `L` must be a single integer between 2 and"
  )
})

test_that("placebo_effects errors for multiple treated units", {
  expect_error(
    placebo_effects(fitN_int, type = "donor"),
    "currently supported only for models with a single treated unit"
  )
})

test_that("placebo_effects errors when no saved data", {
  fit_nosave <- fit1_int
  fit_nosave$data <- NULL
  expect_error(
    placebo_effects(fit_nosave, type = "donor"),
    "does not contain the original data"
  )
})
