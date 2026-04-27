test_that("logistic_normal() and dirichlet() return correct objects", {
  x <- logistic_normal()
  expect_s3_class(x, "omega_prior")
  expect_equal(x$distribution, "logistic_normal")
  expect_null(x$kappa)
  expect_null(x$r_ess)

  x <- logistic_normal(kappa = 2)
  expect_equal(x$kappa, 2)
  expect_null(x$r_ess)

  x <- logistic_normal(r_ess = 0.3)
  expect_null(x$kappa)
  expect_equal(x$r_ess, 0.3)

  x <- dirichlet(kappa = 0.5)
  expect_s3_class(x, "omega_prior")
  expect_equal(x$distribution, "dirichlet")
  expect_equal(x$kappa, 0.5)
})

test_that("omega_prior constructors reject invalid arguments", {
  expect_error(logistic_normal(kappa = -1), "kappa")
  expect_error(logistic_normal(r_ess = 1.5), "r_ess")
  expect_error(logistic_normal(kappa = 1, r_ess = 0.3), "Only one")
  expect_error(dirichlet(kappa = 0, r_ess = 0.3), "Only one")
})

test_that("format.omega_prior produces correct strings", {
  expect_equal(format(logistic_normal(kappa = 2)), "logistic_normal(kappa = 2)")
  expect_equal(format(logistic_normal(r_ess = 0.1)), "logistic_normal(r_ess = 0.1)")
  expect_equal(format(logistic_normal()), "logistic_normal()")
  expect_equal(format(dirichlet(kappa = 0.5)), "dirichlet(kappa = 0.5)")
})

test_that("resolve_kappa() resolves kappa correctly", {
  J <- 10L
  expect_equal(resolve_kappa(logistic_normal(kappa = 3), J), 3)
  expect_equal(
    resolve_kappa(logistic_normal(r_ess = 0.25), J),
    select_kappa(J, 0.25)
  )
  expect_equal(
    resolve_kappa(logistic_normal(), J),
    select_kappa(J, 0.25)
  )
  expect_equal(resolve_kappa(dirichlet(kappa = 1), J), 1)
  expect_error(resolve_kappa(dirichlet(r_ess = 0.3), J), "not yet implemented")
  expect_error(resolve_kappa(dirichlet(), J), "not yet implemented")
})

test_that("bscm() rejects invalid omega_prior", {
  expect_error(
    bscm(y ~ 1, single_treated, "treatment", omega_prior = "logistic_normal"),
    "omega_prior"
  )
})
