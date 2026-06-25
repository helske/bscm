test_that("prior constructors set correct structure", {
  pr <- normal_pr(location = c(0, 1, 3), scale = c(1, 2, 2))
  expect_s3_class(pr, "bscm_prior")
  expect_equal(pr$distribution, "normal")
  expect_equal(pr$length, 3)
  expect_equal(pr$npar, 2)
})
test_that("length mismatch throws error in prior constructors", {
  expect_error(
    normal_pr(location = 0:1, scale = 1:3),
    "match"
  )
})
test_that("scalars are recycled in prior constructors", {
  expect_equal(
    beta_pr(shape1 = 5, shape2 = c(2, 3, 1)),
    structure(
      list(
        distribution = "beta",
        shape1 = rep(5, 3),
        shape2 = c(2, 3, 1),
        npar = 2,
        length = 3L
      ),
      class = "bscm_prior"
    )
  )
  expect_equal(
    gamma_pr(c(2, 3, 1), 0.1),
    structure(
      list(
        distribution = "gamma",
        shape = c(2, 3, 1),
        rate = rep(0.1, 3),
        npar = 2,
        length = 3L
      ),
      class = "bscm_prior"
    )
  )
})
test_that("as.character for prior formats priors correctly", {
  pr <- gamma_pr(2, 1)
  out <- as.character(pr)
  expect_equal(out, "gamma(2, 1)")
  pr <- gamma_pr(1:2, 1:2)
  out <- as.character(pr)
  expect_equal(out, c("gamma(1, 1)", "gamma(2, 2)"))
})
test_that("print method for prior returns invisible object", {
  pr <- exponential_pr(rate = 1)
  expect_invisible(print(pr))
})

test_that("logistic_normal() and dirichlet() return correct objects", {
  x <- logistic_normal_pr(2)
  expect_s3_class(x, "bscm_prior")
  expect_equal(x$distribution, "logistic_normal")
  expect_equal(x$scale, 2)
  z <- 0.5
  x <- dirichlet_pr(z)
  expect_s3_class(x, "bscm_prior")
  expect_equal(x$distribution, "dirichlet")
  expect_equal(x$concentration, z)
})
