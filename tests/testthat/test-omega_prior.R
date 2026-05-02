test_that("logistic_normal() and dirichlet() return correct objects", {
  x <- logistic_normal(2)
  expect_s3_class(x, "omega_prior")
  expect_equal(x$distribution, "logistic_normal")
  expect_equal(x$kappa, 2)
  x <- logistic_normal(kappa = 2)
  expect_equal(x$kappa, 2)
  z <- 0.5
  x <- dirichlet(z)
  expect_s3_class(x, "omega_prior")
  expect_equal(x$distribution, "dirichlet")
  expect_equal(x$kappa, 0.5)
})

test_that("omega_prior constructors reject invalid arguments", {
  expect_error(logistic_normal(kappa = -1), "kappa")
  expect_error(logistic_normal(), "kappa")
})

test_that("bscm() rejects invalid omega_prior", {
  expect_error(
    bscm(y ~ 1, single_treated, "treatment", omega_prior = "logistic_normal"),
    "omega_prior"
  )
})
