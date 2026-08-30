test_that("treatment_effect returns data.frame with correct structure", {
  d <- treatment_effect(fit1_noint)
  expect_s3_class(d, "tbl_df")
  expect_identical(names(d)[1:4], c("id", "time", "treatment", "mean"))
  expect_equal(nrow(d), get_T_total(fit1_noint))

  d <- treatment_effect(fitN_int, average = TRUE)
  expect_equal(d$time, get_times(fitN_int))
  expect_false("id" %in% names(d))

  d <- treatment_effect(fitN_int)
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

test_that("treatment_effect returns unit-specific effects by default", {
  d <- treatment_effect(fitN_int)
  expect_true("id" %in% names(d))
  expect_equal(nrow(d), get_T_total(fitN_int) * get_N(fitN_int))
})

test_that("average_treatment_effect returns correct structure", {
  d <- average_treatment_effect(fit1_noint)
  expect_s3_class(d, "tbl_df")
  expect_identical(names(d)[1:3], c("id", "treatment", "mean"))
  expect_equal(nrow(d), 2L)

  d <- average_treatment_effect(fitN_int, average = FALSE)
  expect_equal(
    nrow(d),
    get_N(fitN_int) * 2L
  )
  expect_true("id" %in% names(d))
})

test_that("average_treatment_effect average argument is ignored for N=1", {
  d1 <- average_treatment_effect(fit1_int)
  d2 <- average_treatment_effect(fit1_int, average = TRUE)
  expect_equal(d1, d2)
})
