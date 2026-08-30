test_that("covariate_imbalance has expected columns", {
  suppressWarnings(d <- covariate_imbalance(fit1_x))
  expect_s3_class(d, "tbl_df")
  expect_named(
    d,
    c(
      "id",
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
  suppressWarnings(d <- covariate_imbalance(fitN_xz))
  expect_s3_class(d, "tbl_df")
  expect_setequal(unique(d$id), get_treated(fitN_xz))
  expect_equal(nrow(d), 120L)
})

test_that("covariate_imbalance returns draws when summary = FALSE", {
  d <- covariate_imbalance(fit1_x, summary = FALSE)
  expect_named(d, c("id", "time", "imbalance"))
  expect_true(posterior::is_rvar(d$imbalance))
})

test_that("covariate_imbalance accepts any number of quantiles", {
  suppressWarnings(d <- covariate_imbalance(fit1_x, probs = numeric(0)))
  expect_false(any(grepl("^q", names(d))))
  suppressWarnings(d <- covariate_imbalance(fit1_x, probs = c(0.05, 0.5, 0.95)))
  expect_true(all(c("q5", "q50", "q95") %in% names(d)))
})

test_that("covariate_imbalance errors for model without predictors", {
  expect_error(
    covariate_imbalance(fit1_int),
    "The model does not contain any predictors"
  )
})

test_that("covariate_imbalance validates its arguments", {
  expect_error(
    covariate_imbalance(fit1_x, probs = "x"),
    "Argument `probs` must be a <numeric> vector"
  )
  expect_error(
    covariate_imbalance(fit1_x, summary = "yes"),
    "Argument `summary` must be a single <logical> value"
  )
})

test_that("plot_covariate_imbalance returns a ggplot per treated unit", {
  expect_s3_class(plot_covariate_imbalance(fit1_x), "ggplot")
  suppressWarnings(p <- plot_covariate_imbalance(fitN_xz))
  expect_named(p, get_treated(fitN_xz))
  suppressWarnings(
    expect_s3_class(plot_covariate_imbalance(fitN_xz, average = TRUE), "ggplot")
  )
  suppressWarnings(
    expect_s3_class(
      plot_covariate_imbalance(fitN_xz, unit = get_treated(fitN_xz)[1]),
      "ggplot"
    )
  )
})

test_that("plot_covariate_imbalance requires two probabilities", {
  expect_error(
    plot_covariate_imbalance(fit1_x, probs = 0.5),
    "Argument `probs` must be a <numeric> vector of length 2"
  )
  expect_error(
    plot_covariate_imbalance(fit1_x, probs = numeric(0)),
    "Argument `probs` must be a <numeric> vector of length 2"
  )
})
