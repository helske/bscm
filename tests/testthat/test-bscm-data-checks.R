test_that("bscm() rejects unbalanced panel data", {
  d <- single_treated[-1L, ]
  expect_error(
    bscm(y ~ 1, data = d, treatment = "treatment"),
    "Data is not balanced"
  )
  d <- single_treated
  d <- d[d$time != 5 | d$id < 30, ]
  d <- d[d$time != 10 | d$id > 30, ]
  expect_error(
    fit <- bscm(
      y ~ x, data = d, 
      treatment = "treatment", algorithm = "Fixed_param", 
      chains = 1, refresh = 0, iter = 2
    ),
    "Data is not balanced"
  )
})

test_that("bscm() rejects gaps in treatment variable", {
  d <- single_treated
  d$treatment[d$id == 1 & d$time == 5] <- 0L
  expect_error(
    bscm(y ~ 1, data = d, treatment = "treatment"),
    "There should be no gaps in the treatment"
  )
})

test_that("bscm() accepts multiple treated units", {
  expect_error(
    bscm(
      y ~ 1, data = multiple_treated, treatment = "treatment", 
      algorithm = "Fixed_param", chains = 1, refresh = 0, iter = 2
    ),
    NA
  )
  d <- multiple_treated
  d$treatment[d$id == 6 & d$time > -2] <- 1
  d$treatment[d$id == 1] <- 0
  expect_error(
    fit <- bscm(
      y ~ 1, data = d, 
      treatment = "treatment", algorithm = "Fixed_param", 
      chains = 1, refresh = 0, iter = 2
    ),
    NA
  )
  expect_equal(
    get_donors(fit), 
    as.character(setdiff(seq_len(max(d$id)), c(2:3, 6)))
  )
})
test_that("bscm() rejects missing values", {
  d <- single_treated
  d$y[10] <- NA
  expect_error(
    bscm(y ~ 1, data = d, treatment = "treatment"),
    "Missing values are not supported\\."
  )
  d <- single_treated
  d$x[10] <- NA
  expect_error(
    bscm(y ~ x, data = d, treatment = "treatment"),
    "Missing values are not supported\\."
  )
  d <- single_treated
  d$time[10] <- NA
  expect_error(
    bscm(y ~ x, data = d, treatment = "treatment"),
    "Time index variable `time` must be of type"
  )
})
test_that("bscm() accepts gaps in time", {
  d <- single_treated
  d <- d[d$time != 5, ]
  expect_error(
    fit <- bscm(
      y ~ x, data = d, 
      treatment = "treatment", algorithm = "Fixed_param", 
      chains = 1, refresh = 0, iter = 2
    ),
    NA
  )
})
test_that("bscm() rejects constant outcome in pre-treatment", {
  d <- single_treated
  d$y[d$id == 1 & d$treatment == 0] <- 5
  expect_error(
    bscm(y ~ 1, data = d, treatment = "treatment"),
    "Outcome variable cannot be constant in the pre-treatment period"
  )
})

test_that("bscm() rejects constant predictor for all units", {
  d <- single_treated
  d$x <- 1
  expect_warning(
    bscm(y ~ x, data = d, treatment = "treatment", 
         algorithm = "Fixed_param", chains = 1, iter = 1, refresh = 0),
    "Model has unit-specific intercepts and predictors which do not vary"
  )
})

test_that("bscm() validates effective_donors range", {
  expect_error(
    bscm(
      y ~ 1, data = single_treated, treatment = "treatment",
      effective_donors = 30L
    ),
    "Argument `effective_donors` should be between 2 and 29\\."
  )
})
