test_that("as_draws returns a draws_array", {
  d <- as_draws(fit1_int, "alpha")
  expect_s3_class(d, "draws_array")
  expect_equal(dim(d), c(ndraws(fit1_int) / 2L, 2L, 1L))
  d <- as_draws(fit1_int, c("alpha", "y_mean"))
  expect_s3_class(d, "draws_array")
  expect_equal(dim(d), c(ndraws(fit1_int) / 2L, 2L, 41L))
  expect_equal(
    posterior::variables(d),
    c("alpha[1]", paste0("y_mean[", seq_len(40), ",1]"))
  )
  d <- as_draws(
    fit1_int, 
    c("omega", "effect", "y_mean", "y_rep", "log_lik"), 
    include = FALSE
  )
  expect_equal(
    posterior::variables(d),
    c(
      "sigma[1]", "alpha[1]", "R2[1]", 
      "RMSE_pre[1]", "RMSE_post[1]", "RMSE_ratio[1]", "effective_donors[1]", 
      "avg_effect_pre", "avg_effect_post", "avg_RMSE_pre", "avg_RMSE_post",
      "avg_RMSE_ratio", "avg_effective_donors", "lp__"
    )
  )
})
