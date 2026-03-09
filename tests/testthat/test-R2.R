test_that("bayes_R2 returns a numeric vector", {
  r2 <- bayes_R2(fit_intercept, summary = FALSE)
  expect_type(r2, "double")
  expect_equal(length(r2), ndraws(fit_intercept))
  expect_true(all(r2 >= -1 & r2 <= 1))
})

test_that("bayes_R2 returns a summary", {
  r2 <- bayes_R2(fit_intercept, summary = TRUE)
  expect_s3_class(r2, "draws_summary")
})

test_that("loo_R2 returns a numeric vector", {
  suppressWarnings(r2 <- loo_R2(fit_x, summary = FALSE))
  expect_type(r2, "double")
  expect_equal(length(r2), ndraws(fit_x))
  expect_true(all(r2 >= -1 & r2 <= 1))
})
