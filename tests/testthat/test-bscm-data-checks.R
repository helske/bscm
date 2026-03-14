test_that("bscm() rejects unbalanced panel data", {
  d <- simulated_data[-1L, ]
  expect_error(
    bscm(y ~ 1, data = d, treatment = "treatment"),
    "Data is not balanced\\. All units should have equal number of time points\\."
  )
})

test_that("bscm() rejects multiple treated units", {
  d <- simulated_data
  d$treatment[d$id == 2 & d$time >= 0] <- 1L
  expect_error(
    bscm(y ~ 1, data = d, treatment = "treatment"),
    "Only the case of a single treated unit is currently implemented\\."
  )
})

test_that("bscm() rejects gaps in treatment variable", {
  d <- simulated_data
  d$treatment[d$id == 1 & d$time == 5] <- 0L
  expect_error(
    bscm(y ~ 1, data = d, treatment = "treatment"),
    "There should be no gaps in the treatment\\. Check the treatment variable\\."
  )
})

test_that("bscm() rejects missing values in outcome", {
  d <- simulated_data
  d$y[1] <- NA
  expect_error(
    bscm(y ~ 1, data = d, treatment = "treatment"),
    "Missing values are not supported\\."
  )
})

test_that("bscm() rejects constant outcome in pre-treatment", {
  d <- simulated_data
  d$y[d$id == 1 & d$treatment == 0] <- 5
  expect_error(
    bscm(y ~ 1, data = d, treatment = "treatment"),
    "Outcome variable cannot be constant in the pre-treatment period\\."
  )
})

test_that("bscm() rejects constant predictor for all units", {
  d <- simulated_data
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
      y ~ 1, data = simulated_data, treatment = "treatment",
      effective_donors = 30L
    ),
    "Argument `effective_donors` should be between 2 and 29\\."
  )
})

test_that("bscm() warns about small kappa from effective_donors", {
  expect_warning(
      bscm(
        y ~ 1, data = simulated_data, treatment = "treatment",
        effective_donors = 2L, chains = 1, algorithm = "Fixed_param",
        refresh = 0, mcmc_diagnostics = FALSE
      ),
    "implies a Dirichlet prior"
  )
})
