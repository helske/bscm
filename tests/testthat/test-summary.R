test_that("summary returns posterior summaries of the model parameters", {
  suppressWarnings(s <- summary(fit1_int))
  expect_s3_class(s, "tbl_df")
  expect_true(
    all(
      c("variable", "mean", "sd", "q2.5", "q97.5", "rhat") %in% names(s)
    )
  )
  expect_setequal(s$variable, c("alpha[1]", "sigma[1]"))
})

test_that("summary respects the parameters argument", {
  suppressWarnings(s <- summary(fit1_x, parameters = "beta"))
  expect_equal(s$variable, "beta[1]")
})

test_that("summary validates probs argument", {
  expect_error(
    summary(fit1_noint, probs = "x"),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("summary respects custom probs", {
  suppressWarnings(s <- summary(fit1_int, probs = c(0.1, 0.9)))
  expect_true(all(c("q10", "q90") %in% names(s)))
})
