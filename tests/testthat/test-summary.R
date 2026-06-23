test_that("summary returns a bscmfit_summary", {
  suppressWarnings(s <- summary(fit1_int))
  expect_s3_class(s, "summary_bscmfit")
  expect_equal(
    s$variable, 
    c(
      "Intercept", 
      "Residual SD",
      "Bayesian R2",
      "Effective number of donors",
      "Average pre-treatment effect", 
      "Average post-treatment effect", 
      "Pre-treatment RMSE", 
      "Post-treatment RMSE"
    )
  )
})

test_that("summary validates probs argument", {
  expect_error(
    summary(fit1_noint, probs = "x"),
    "Argument `probs` must be a <numeric> vector"
  )
})
