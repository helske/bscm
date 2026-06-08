test_that("effective_donors returns data.frame", {
  d <- effective_donors(fit1_int, probs = 0.3)
  expect_s3_class(d, "tbl_df")
  expect_named(
    d,
    c(
      "mean",
      "sd",
      "q30",
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean"
    )
  )
  expect_equal(nrow(d), 1L)

  d <- effective_donors(fitN_int)
  expect_s3_class(d, "tbl_df")
  expect_named(
    d,
    c(
      "id",
      "mean",
      "sd",
      "q2.5",
      "q97.5",
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean"
    )
  )
  d <- effective_donors(fitN_int, average = TRUE)
  expect_s3_class(d, "tbl_df")
  expect_null(d[["variable"]])
})

test_that("effective_donors respects custom probs", {
  d <- effective_donors(fit1_int, probs = c(0.1, 0.9, 0.99))
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
  expect_true("q99" %in% names(d))
  d <- effective_donors(fit1_int, probs = 0.5)
  expect_true("q50" %in% names(d))
  d <- effective_donors(fit1_int, probs = numeric(0))
  expect_named(
    d,
    c(
      "mean",
      "sd",
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean"
    )
  )
})

test_that("effective_donors validates probs argument", {
  expect_error(
    effective_donors(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("effective_donors average argument works", {
  expect_error(
    effective_donors(fit1_int, average = "yes"),
    "Argument `average` must be a single"
  )
  d1 <- effective_donors(fit1_int)
  d2 <- effective_donors(fit1_int, average = FALSE)
  expect_equal(d1, d2)
  d1 <- effective_donors(fitN_int)
  d2 <- effective_donors(fitN_int, average = FALSE)
  expect_equal(d1, d2)
})
