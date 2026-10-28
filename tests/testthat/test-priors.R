test_that("prior constructors set correct structure", {
  pr <- normal_pr(location = c(0, 1, 3), scale = c(1, 2, 2))
  expect_s3_class(pr, "bscm_prior")
  expect_equal(pr$distribution, "normal")
  expect_equal(pr$length, 3)
  expect_equal(pr$npar, 2)
})
test_that("student_pr sets correct structure", {
  pr <- student_pr(df = 4, location = 0, scale = 1)
  expect_s3_class(pr, "bscm_prior")
  expect_equal(pr$distribution, "student_t")
  expect_equal(pr$npar, 3)
  expect_equal(pr$length, 1)
})
test_that("half_normal_pr sets correct structure", {
  pr <- half_normal_pr(scale = 2)
  expect_s3_class(pr, "bscm_prior")
  expect_equal(pr$distribution, "half_normal")
  expect_equal(pr$npar, 1)
  expect_equal(pr$scale, 2)
})
test_that("prior constructors reject non-positive parameters", {
  expect_error(normal_pr(0, -1), "positive")
  expect_error(gamma_pr(-1, 1), "positive")
  expect_error(exponential_pr(-1), "positive")
  expect_error(half_normal_pr(-1), "positive")
  expect_error(beta_pr(-1, 1), "positive")
  expect_error(dirichlet_pr(-1), "positive")
  expect_error(logistic_normal_pr(-1), "positive")
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
test_that("as.character works for all prior families", {
  expect_equal(as.character(normal_pr(0, 1)), "normal(0, 1)")
  expect_equal(as.character(student_pr(4, 0, 1)), "student(4, 0, 1)")
  expect_equal(as.character(beta_pr(2, 3)), "beta(2, 3)")
  expect_equal(as.character(exponential_pr(2)), "exponential(2)")
  expect_equal(as.character(half_normal_pr(2)), "half_normal(2)")
  expect_equal(as.character(dirichlet_pr(1)), "dirichlet(1)")
  expect_equal(as.character(logistic_normal_pr(2)), "logistic_normal(2)")
  expect_equal(as.character(cauchy_pr(0, 1)), "cauchy(0, 1)")
  expect_equal(as.character(laplace_pr(0, 1)), "laplace(0, 1)")
  expect_equal(as.character(inv_gamma_pr(2, 1)), "inv_gamma(2, 1)")
  expect_equal(
    as.character(scaled_inv_chi_pr(4, 1)),
    "scaled_inv_chi(4, 1)"
  )
  expect_equal(
    as.character(sqrt_inv_gamma_pr(2, 1)),
    "sqrt_inv_gamma(2, 1)"
  )
  expect_equal(
    as.character(logistic_laplace_pr(2)),
    "logistic_laplace(2)"
  )
})
test_that("inv_gamma_pr uses shape and scale", {
  pr <- inv_gamma_pr(shape = 2, scale = c(1, 3))
  expect_equal(pr$shape, c(2, 2))
  expect_equal(pr$scale, c(1, 3))
  expect_null(pr$rate)
})
test_that("prior constructors check the correct arguments", {
  expect_error(student_pr(-1, 0, 1), "df")
  expect_error(sqrt_inv_gamma_pr(-1, 1), "shape")
  expect_error(scaled_inv_chi_pr(4, -1), "scale")
})
test_that("available_priors lists all prior constructors", {
  pr <- available_priors()
  expect_s3_class(pr, "data.frame")
  expect_equal(pr$distribution, sort(pr$distribution))
  fns <- sub("()", "", pr$constructor, fixed = TRUE)
  expect_true(all(vapply(fns, exists, logical(1))))
  expect_true(all(pr$support %in% c("real", "positive", "(0, 1)", "simplex")))
})
test_that("prior_to_stan has a code for all available priors", {
  pr <- available_priors()
  for (i in seq_len(nrow(pr))) {
    x <- structure(
      list(
        distribution = pr$distribution[i],
        location = 0,
        scale = 1,
        shape = 1,
        rate = 1,
        df = 1,
        shape1 = 1,
        shape2 = 1,
        concentration = 1,
        length = 1L,
        zeta = FALSE
      ),
      class = "bscm_prior"
    )
    type <- if (pr$support[i] == "simplex") "omega" else "sigma"
    expect_false(anyNA(prior_to_stan(x, type)[[1]]))
  }
})
test_that("get_priors returns a list of bscm_prior objects", {
  pr_list <- get_priors(fit1_int, FALSE)
  expect_type(pr_list, "list")
  expect_true(all(vapply(pr_list, inherits, logical(1L), "bscm_prior")))
})
test_that("print method for prior returns invisible object", {
  pr <- exponential_pr(rate = 1)
  expect_invisible(print(pr))
})

test_that("logistic_normal() and dirichlet() return correct objects", {
  x <- logistic_normal_pr(2)
  expect_s3_class(x, "bscm_prior")
  expect_equal(x$distribution, "logistic_normal")
  expect_equal(x$concentration, 2)
  z <- 0.5
  x <- dirichlet_pr(z)
  expect_s3_class(x, "bscm_prior")
  expect_equal(x$distribution, "dirichlet")
  expect_equal(x$concentration, z)
})
