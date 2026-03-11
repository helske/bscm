test_that("summary returns a bscmfit_summary list", {
  suppressWarnings(s <- summary(fit_intercept))
  expect_type(s, "list")
  expected_elements <- c(
    "effects", "synthetic", "weights", "cumulative_effects", 
    "relative_change", "parameters", "RMSE", "R2", 
    "effective_donors", "average_effects"
  )
  expect_true(all(expected_elements %in% names(s)))
  d <- s$effects
  expect_s3_class(d, "data.frame")
  expect_true("time" %in% names(d))
  expect_true("mean" %in% names(d))
  expect_true("sd" %in% names(d))
  expect_equal(nrow(d), fit_intercept$setup$T_total)
})

test_that("summary with include argument selects subset", {
  suppressWarnings(s <- summary(fit_x, include = "effects"))
  expect_true("effects" %in% names(s))
  expect_false("weights" %in% names(s))
})

test_that("summary validates probs argument", {
  expect_error(
    summary(fit_intercept, probs = "x"),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("summary rejects invalid include values", {
  expect_error(
    summary(fit_intercept, include = "dfsdg"),
    "Argument `include` contains invalid values"
  )
})

