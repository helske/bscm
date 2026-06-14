test_that("synthetic_control returns data.frame with correct structure for single treated", {
  d <- synthetic_control(fit1_int)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), get_T_total(fit1_int))
  expect_true(all(c("time", "treatment", "mean", "sd") %in% names(d)))
})

test_that("synthetic_control returns raw draws when summary = FALSE", {
  d <- synthetic_control(fit1_int, summary = FALSE)
  expect_s3_class(d, "data.frame")
  expect_true(posterior::is_rvar(d$y_rep))
  expect_equal(nrow(d), get_T_total(fit1_int))
})

test_that("synthetic_control respects custom probs", {
  d <- synthetic_control(fit1_int, probs = c(0.1, 0.9))
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
})

test_that("synthetic_control validates summary argument", {
  expect_error(
    synthetic_control(fit1_int, summary = "yes"),
    "Argument `summary` must be a single"
  )
})

test_that("synthetic_control handles N > 1 with unit column", {
  d <- synthetic_control(fitN_int)
  expect_equal(nrow(d), get_T_total(fitN_int) * get_N(fitN_int))
  expect_true("id" %in% names(d))
  expect_setequal(unique(d$id), get_treated(fitN_int))
})

test_that("synthetic_control errors when model has no saved data", {
  fit_nosave <- fit1_int
  fit_nosave$data <- NULL
  expect_error(
    synthetic_control(fit_nosave),
    "does not contain the original data"
  )
})
