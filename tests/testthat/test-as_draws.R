test_that("as_draws returns a draws_array", {
  d <- as_draws(fit_intercept, "alpha")
  expect_s3_class(d, "draws_array")
  expect_equal(dim(d), c(ndraws(fit_intercept), 1L, 1L))
  d <- as_draws(fit_intercept, c("alpha", "synthetic_mean"))
  expect_s3_class(d, "draws_array")
  expect_equal(dim(d), c(ndraws(fit_intercept), 1L, 51L))
  expect_equal(
    posterior::variables(d),
    c("alpha", paste0("synthetic_mean[", seq_len(50),"]"))
  )
  d <- as_draws(
    fit_intercept, 
    c("omega", "effect", "synthetic_mean", "synthetic_y", "log_lik",
      "avg_effect_post_cumulative", "relative_change"), 
    include = FALSE
  )
  expect_equal(
    posterior::variables(d),
    c(
      "sigma", "alpha", "R2", "avg_effect_pre", "avg_effect_post", 
      "RMSE_pre", "RMSE_post", "RMSE_ratio", "effective_donors", "lp__"
    )
  )
})
