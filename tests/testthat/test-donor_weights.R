test_that("donor_weights returns tibble with correct structure for single treated", {
  d <- donor_weights(fit1_int)
  expect_s3_class(d, "tbl_df")
  expect_named(
    d,
    c(
      "treated",
      "donor",
      "mean",
      "sd",
      "q2.5",
      "q50",
      "q97.5",
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean"
    )
  )
  expect_equal(nrow(d), get_J(fit1_int))
})

test_that("donor_weights returns raw draws when summary = FALSE", {
  d <- donor_weights(fit1_int, summary = FALSE)
  expect_s3_class(d, "tbl_df")
  expect_named(d, c("treated", "donor", "weight"))
  expect_true(posterior::is_rvar(d$weight))
  expect_equal(nrow(d), get_J(fit1_int))
})

test_that("donor_weights respects custom probs", {
  d <- donor_weights(fit1_int, probs = c(0.1, 0.9))
  expect_true("q10" %in% names(d))
  expect_true("q90" %in% names(d))
})

test_that("donor_weights validates summary argument", {
  expect_error(
    donor_weights(fit1_int, summary = "yes"),
    "Argument `summary` must be a single"
  )
})

test_that("donor_weights validates probs argument", {
  expect_error(
    donor_weights(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("donor_weights covers all treated units for N > 1", {
  suppressWarnings(d <- donor_weights(fitN_int))
  expect_setequal(unique(d$treated), get_treated(fitN_int))
  expect_equal(nrow(d), get_N(fitN_int) * get_J(fitN_int))
})

test_that("donor_weights sum to 1 per treated unit", {
  d <- donor_weights(fit1_int)
  expect_equal(sum(d$mean), 1, tolerance = 0.01)
})
