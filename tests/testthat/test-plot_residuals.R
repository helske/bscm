test_that("plot_residuals returns ggplot for single treated unit", {
  p <- plot_residuals(fit1_int)
  expect_s3_class(p, "ggplot")
})

test_that("plot_residuals returns named list of ggplots for N > 1", {
  suppressWarnings(plots <- plot_residuals(fitN_int))
  expect_type(plots, "list")
  expect_named(plots, get_treated(fitN_int))
  expect_true(all(vapply(plots, inherits, logical(1L), "ggplot")))
})

test_that("plot_residuals returns a single ggplot for a selected unit", {
  suppressWarnings(
    p <- plot_residuals(fitN_int, unit = get_treated(fitN_int)[1])
  )
  expect_s3_class(p, "ggplot")
  suppressWarnings(
    p <- plot_residuals(
      fitN_int,
      type = "autocorrelation",
      unit = get_treated(fitN_int)[2]
    )
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_residuals type = 'autocorrelation' returns ggplot", {
  p <- plot_residuals(fit1_int, type = "autocorrelation")
  expect_s3_class(p, "ggplot")
})

test_that("plot_residuals validates type argument", {
  expect_error(
    plot_residuals(fit1_int, type = "invalid"),
    "Argument `type` must be either"
  )
})

test_that("plot_residuals validates unit argument", {
  expect_error(
    suppressWarnings(plot_residuals(fitN_int, unit = "not_a_unit")),
    "Argument `unit` must be one of the treated units"
  )
})

test_that("plot_residuals validates probs argument length", {
  expect_error(
    plot_residuals(fit1_int, probs = c(0.025, 0.5, 0.975)),
    "Argument `probs` must be a"
  )
})

test_that("residual_acf returns summaries with correct structure", {
  d <- residual_acf(fit1_int, max_lag = 4L)
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 4L)
  expect_equal(d$lag, seq_len(4L))
  expect_named(
    d,
    c(
      "id",
      "lag",
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
})

test_that("residual_acf returns draws when summary = FALSE", {
  d <- residual_acf(fit1_int, summary = FALSE)
  expect_named(d, c("id", "lag", "ac"))
  expect_true(posterior::is_rvar(d$ac))
})

test_that("residual_acf is computed separately for each treated unit", {
  suppressWarnings(d <- residual_acf(fitN_int, max_lag = 3L))
  expect_setequal(unique(d$id), get_treated(fitN_int))
  expect_equal(nrow(d), 3L * get_N(fitN_int))
})

test_that("residual_acf accepts any number of quantiles", {
  d <- residual_acf(fit1_int, probs = numeric(0))
  expect_false(any(grepl("^q", names(d))))
  d <- residual_acf(fit1_int, probs = c(0.05, 0.5, 0.95))
  expect_true(all(c("q5", "q50", "q95") %in% names(d)))
})

test_that("residual_acf validates max_lag argument", {
  expect_error(
    residual_acf(fit1_int, max_lag = -1L),
    "Argument `max_lag` must be a single positive integer"
  )
  expect_error(
    residual_acf(fit1_int, max_lag = 1.5),
    "Argument `max_lag` must be a single positive integer"
  )
})
