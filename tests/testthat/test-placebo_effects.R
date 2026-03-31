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
