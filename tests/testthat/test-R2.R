test_that("bayes_R2 returns draws_summary", {
  r2 <- bayes_R2(fit1_int, probs = c(0.05, 0.5, 0.9, 0.99))
  expect_s3_class(r2, "draws_summary")
  expect_identical(dim(r2), c(1L, 11L))

  r2 <- bayes_R2(fitN_xz, probs = numeric(0))
  expect_s3_class(r2, "draws_summary")
  expect_identical(dim(r2), c(3L, 7L))
  expect_named(
    r2,
    c("id", "variable", "mean", "sd", "rhat", "ess_bulk", "ess_tail")
  )
})

test_that("loo_R2 returns draws_summary", {
  suppressWarnings(r2 <- loo_R2(fit1_x, probs = 0.2))
  expect_s3_class(r2, "draws_summary")
  expect_identical(dim(r2), c(1L, 8L))
})
