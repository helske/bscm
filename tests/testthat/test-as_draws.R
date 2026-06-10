test_that("as_draws returns a draws_array", {
  d <- as_draws(fit1_int, "alpha")
  expect_s3_class(d, "draws_array")
  expect_equal(dim(d), c(ndraws(fit1_int) / 2L, 2L, 1L))
  d <- as_draws(fit1_int, c("alpha", "y_mean"))
  expect_s3_class(d, "draws_array")
  T_ <- get_T_total(fit1_int)
  expect_equal(dim(d), c(ndraws(fit1_int) / 2L, 2L, 1 + T_))
  expect_equal(
    posterior::variables(d),
    c("alpha[1]", paste0("y_mean[", seq_len(T_), ",1]"))
  )
  d <- as_draws(
    fit1_int,
    c("omega", "y_mean", "y_rep"),
    include = FALSE
  )
  expect_equal(
    posterior::variables(d),
    c("sigma[1]", "alpha[1]", "lp__")
  )
})
