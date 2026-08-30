test_that("covariate_adjustment has expected columns", {
  suppressWarnings(d <- covariate_adjustment(fit1_x))
  expect_s3_class(d, "tbl_df")
  expect_named(
    d,
    c(
      "id",
      "variable",
      "time",
      "mean",
      "sd",
      "q2.5",
      "q97.5",
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean"
    )
  )
  expect_setequal(unique(d$variable), fit1_x$setup$beta_names)
})

test_that("covariate_adjustment returns draws when summary = FALSE", {
  d <- covariate_adjustment(fit1_x, summary = FALSE)
  expect_named(d, c("id", "variable", "time", "adjustment"))
  expect_true(posterior::is_rvar(d$adjustment))
})

test_that("covariate_adjustment accepts any number of quantiles", {
  suppressWarnings(d <- covariate_adjustment(fit1_x, probs = numeric(0)))
  expect_false(any(grepl("^q", names(d))))
  suppressWarnings(
    d <- covariate_adjustment(fit1_x, probs = c(0.05, 0.5, 0.95))
  )
  expect_true(all(c("q5", "q50", "q95") %in% names(d)))
})

test_that("covariate_adjustment is unit specific unless averaged", {
  suppressWarnings(d <- covariate_adjustment(fitN_tv))
  expect_setequal(unique(d$id), get_treated(fitN_tv))
  suppressWarnings(d_avg <- covariate_adjustment(fitN_tv, average = TRUE))
  expect_false("id" %in% names(d_avg))
  expect_equal(nrow(d_avg) * get_N(fitN_tv), nrow(d))
})

test_that("covariate_adjustment handles fixed and time-varying coefficients", {
  suppressWarnings(d <- covariate_adjustment(fitN_tv))
  expect_setequal(unique(d$variable), fitN_tv$setup$beta_names)
  # predictors with a fixed coefficient must not be dropped to NA in a model
  # that also contains time-varying coefficients
  expect_false(anyNA(d$mean))
})

test_that("covariate_adjustment errors for model without predictors", {
  expect_error(
    covariate_adjustment(fit1_int),
    "The model does not contain any predictors"
  )
})

test_that("plot_covariate_adjustment returns a ggplot", {
  expect_s3_class(plot_covariate_adjustment(fit1_x), "ggplot")
  expect_error(
    plot_covariate_adjustment(fit1_x, probs = 0.5),
    "Argument `probs` must be a <numeric> vector of length 2"
  )
  expect_error(
    plot_covariate_adjustment(fit1_x, alpha = 2),
    "Argument `alpha` must be a single <numeric> value"
  )
})
