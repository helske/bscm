test_that("lfo validates L argument is too small", {
  expect_error(
    lfo(fit1_int, L = 1L),
    "Argument `L` must be a single integer between 2 and"
  )
})

test_that("lfo validates L argument is too large", {
  T_pre_min <- min(get_T_pre(fit1_int))
  expect_error(
    lfo(fit1_int, L = T_pre_min - 1L),
    "Argument `L` must be a single integer between 2 and"
  )
})

test_that("lfo validates L argument is an integer", {
  expect_error(
    lfo(fit1_int, L = "a"),
    "Argument `L` must be a single integer between 2 and"
  )
})
