test_that("treatment_effect returns data.frame with correct structure", {
  d <- treatment_effect(fit1_noint)
  expect_s3_class(d, "data.frame")
  expect_identical(names(d)[1:4], c("id", "time", "variable", "mean"))
  expect_equal(nrow(d), get_T_total(fit1_noint))

  d <- treatment_effect(fitN_int)
  expect_equal(d$time, get_times(fitN_int))
  expect_identical(d$variable[1], "Average treatment effect")

  d <- treatment_effect(fitN_int, average = FALSE)
  expect_equal(d$time, rep(get_times(fitN_int), 3))
  expect_equal(d$id, rep(get_treated(fitN_int), each = get_T_total(fitN_int)))
})

test_that("treatment_effect respects custom probs", {
  d <- treatment_effect(fit1_int, probs = c(0.1, 0.9))
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
})

test_that("treatment_effect validates probs argument", {
  expect_error(
    treatment_effect(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("treatment_effect average argument is ignored for N=1", {
  d1 <- treatment_effect(fit1_int)
  d2 <- treatment_effect(fit1_int, average = TRUE)
  expect_equal(d1, d2)
})
