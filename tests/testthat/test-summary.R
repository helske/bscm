test_that("summary returns posterior summaries of the model parameters", {
  suppressWarnings(s <- summary(fit1_int))
  expect_s3_class(s, "tbl_df")
  expect_named(
    s,
    c(
      "parameter", "term", "id", "mean", "sd", "q2.5", "q97.5", "rhat",
      "ess_bulk", "ess_tail", "mcse_mean"
    )
  )
  expect_equal(
    s$parameter,
    c("intercept", "sigma", "concentration", rep("tilt", 3))
  )
  expect_equal(s$term, c(rep(NA, 3), "PC1", "PC2", "PC3"))
  expect_equal(s$id, rep(get_treated(fit1_int), 6))
})

test_that("summary respects the parameters argument", {
  suppressWarnings(s <- summary(fit1_x, parameters = "coef"))
  expect_equal(s$parameter, "coef")
  expect_equal(s$term, "x")
  expect_null(s[["id"]])
  suppressWarnings(s <- summary(fitN_tv, parameters = c("tv_sd", "coef")))
  expect_equal(s$parameter, c("coef", "coef", "tv_sd"))
  expect_equal(s$term, c("z", "x", "x"))
})

test_that("summary returns the unit-specific parameters of each unit", {
  suppressWarnings(s <- summary(fitN_xz, parameters = "sigma"))
  expect_equal(s$id, get_treated(fitN_xz))
  expect_null(s[["term"]])
})

test_that("summary validates the parameters argument", {
  expect_error(
    summary(fit1_int, parameters = "beta"),
    "Argument `parameters` must contain parameters of the model"
  )
  expect_error(
    summary(fit1_int, parameters = "weights"),
    "Argument `parameters` must contain parameters of the model"
  )
})

test_that("summary shows the tilt coefficients by direction", {
  fit <- suppressWarnings(
    bscm(
      y ~ 1,
      data = single_treated,
      priors = list(
        weights = logistic_normal_pr(2, tilt = list(R = 2, scale = 1))
      ),
      chains = 1,
      iter = 200,
      refresh = 0,
      seed = 1,
      mcmc_diagnostics = FALSE
    )
  )
  s <- summary(fit, parameters = "tilt", diagnostics = FALSE)
  expect_equal(s$term, c("PC1", "PC2"))
  expect_equal(s$id, rep(get_treated(fit), 2))
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
