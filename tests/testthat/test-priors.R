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
  expect_equal(as.character(dirichlet_pr(1)), "dirichlet(concentration = 1)")
  expect_equal(
    as.character(logistic_normal_pr(2)),
    "logistic_normal(concentration = 2)"
  )
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
    "logistic_laplace(concentration = 2)"
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
        length = 1L
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
  expect_type(x$concentration, "double")
  expect_null(x$tilt)
})

test_that("simplex priors define fixed and estimated concentrations", {
  x <- logistic_normal_pr(c(1, 2))
  expect_type(x$concentration, "double")
  expect_equal(x$length, 2L)
  x <- logistic_laplace_pr()
  expect_null(x$concentration)
  expect_equal(x$length, 1L)
  x <- dirichlet_pr(inv_gamma_pr(c(2, 3), 1))
  expect_s3_class(x$concentration, "bscm_prior")
  expect_equal(x$concentration, inv_gamma_pr(c(2, 3), 1))
  expect_equal(x$length, 2L)
})

test_that("simplex priors reject invalid concentrations", {
  expect_error(
    dirichlet_pr(NA),
    "Argument `concentration` must be positive number\\(s\\)"
  )
  expect_error(logistic_normal_pr(0:1), "positive")
  expect_error(logistic_normal_pr("1"), "positive")
})

test_that("simplex priors store the tilt", {
  x <- logistic_normal_pr(tilt = list(scale = 2, R = 3))
  expect_equal(x$tilt, list(R = 3L, V = NULL, scale = 2))
  V <- matrix(rnorm(6), 3, dimnames = list(c("a", "b", "c"), NULL))
  x <- dirichlet_pr(1, tilt = list(V = V, scale = 1))
  expect_equal(x$tilt, list(R = 2L, V = V, scale = 1))
})

test_that("simplex priors reject invalid tilts", {
  V <- matrix(1:6, 3, dimnames = list(c("a", "b", "c"), NULL))
  expect_error(
    logistic_normal_pr(tilt = 3),
    "Argument `tilt` must be a named list or `NULL`\\."
  )
  expect_error(
    logistic_normal_pr(tilt = list(R = 2, scale = 1, K = 1)),
    "Invalid element in `tilt`: \"K\"\\."
  )
  expect_error(
    logistic_normal_pr(tilt = list(R = 2, V = V, scale = 1)),
    "must contain either element \"R\" or element \"V\""
  )
  expect_error(
    logistic_normal_pr(tilt = list(scale = 1)),
    "must contain either element \"R\" or element \"V\""
  )
  expect_error(
    logistic_normal_pr(tilt = list(R = 2)),
    "Element \"scale\" of `tilt` must be a single positive number\\."
  )
  expect_error(
    logistic_normal_pr(tilt = list(R = 2, scale = 0)),
    "Element \"scale\" of `tilt` must be a single positive number\\."
  )
  expect_error(
    logistic_normal_pr(tilt = list(R = 1.5, scale = 1)),
    "Element \"R\" of `tilt` must be a single positive integer\\."
  )
  expect_error(
    logistic_normal_pr(tilt = list(R = 0, scale = 1)),
    "Element \"R\" of `tilt` must be a single positive integer\\."
  )
  expect_error(
    logistic_normal_pr(tilt = list(V = unname(V), scale = 1)),
    "Element \"V\" of `tilt` must have unique row names\\."
  )
  expect_error(
    logistic_normal_pr(tilt = list(V = V[c(1, 1, 2), ], scale = 1)),
    "Element \"V\" of `tilt` must have unique row names\\."
  )
  V[1, 1] <- NA
  expect_error(
    logistic_normal_pr(tilt = list(V = V, scale = 1)),
    "Element \"V\" of `tilt` must be a numeric matrix with finite"
  )
})

test_that("as.character shows the concentration and the tilt", {
  expect_equal(
    as.character(logistic_normal_pr()),
    "logistic_normal(concentration = NULL)"
  )
  expect_equal(
    as.character(dirichlet_pr(gamma_pr(2, 1))),
    "dirichlet(concentration = gamma(2, 1))"
  )
  expect_equal(
    as.character(logistic_laplace_pr(tilt = list(R = 3, scale = 2))),
    "logistic_laplace(concentration = NULL, tilt = list(R = 3, scale = 2))"
  )
  V <- matrix(1:6, 3, dimnames = list(c("a", "b", "c"), NULL))
  expect_equal(
    as.character(logistic_normal_pr(2, tilt = list(V = V, scale = 1))),
    "logistic_normal(concentration = 2, tilt = list(R = 2, scale = 1))"
  )
})

test_that("define_priors resolves the prior of the concentrations", {
  setup <- list(
    treated = c("t1", "t2"),
    donors = letters[1:5],
    has_icpt = TRUE,
    has_x = FALSE,
    has_w = FALSE,
    has_ar1 = FALSE,
    marginalize_sigma = TRUE
  )
  d <- list(sd_e = c(1, 2), mean_y = c(0, 1))
  pr <- define_priors(NULL, d, setup, NULL)
  expect_named(pr, c("intercept", "sigma", "weights"))
  expect_s3_class(pr$weights$concentration, "bscm_prior")
  expect_equal(pr$weights$length, 2L)
  expect_equal(pr$weights$concentration, default_zeta_prior(2, 5))
  expect_equal(pr$weights$tilt, list(R = 3L, V = NULL, scale = 2))
  # the default prior has no tilt without marginalized sigma
  setup$marginalize_sigma <- FALSE
  pr <- define_priors(NULL, d, setup, NULL)
  expect_null(pr$weights$tilt)
  setup$marginalize_sigma <- TRUE
  # default prior of the concentration is filled in for a user-defined prior
  tilt <- list(R = 2, scale = 1)
  pr <- define_priors(
    list(weights = dirichlet_pr(tilt = tilt)),
    d,
    setup,
    NULL
  )
  expect_equal(pr$weights$distribution, "dirichlet")
  expect_equal(pr$weights$concentration, default_zeta_prior(2, 5))
  expect_equal(pr$weights$tilt, list(R = 2L, V = NULL, scale = 1))
  # scalar priors are recycled to the number of treated units
  pr <- define_priors(
    list(weights = logistic_normal_pr(gamma_pr(2, 1))),
    d,
    setup,
    NULL
  )
  expect_equal(pr$weights$concentration, gamma_pr(c(2, 2), 1))
  expect_equal(pr$weights$length, 2L)
  pr <- define_priors(list(weights = logistic_normal_pr(3)), d, setup, NULL)
  expect_type(pr$weights$concentration, "double")
  expect_equal(pr$weights$concentration, c(3, 3))
  expect_equal(pr$weights$length, 2L)
  expect_error(
    define_priors(
      list(weights = logistic_normal_pr(gamma_pr(1:3, 1))),
      d,
      setup,
      NULL
    ),
    "Length mismatch for parameter 'concentration'\\."
  )
  expect_error(
    define_priors(list(weights = logistic_normal_pr(1:3)), d, setup, NULL),
    "Length mismatch for parameter 'weights'\\."
  )
  expect_error(
    define_priors(
      list(weights = logistic_normal_pr(beta_pr(2, 2))),
      d,
      setup,
      NULL
    ),
    "Invalid prior family 'beta' for the concentration"
  )
  # the names of the parameters in the Stan model are not accepted
  for (nm in c("omega", "alpha", "beta", "kappa", "rho")) {
    pr <- stats::setNames(list(normal_pr(0, 1)), nm)
    expect_error(
      define_priors(pr, d, setup, NULL),
      paste0("Invalid parameter name in `priors`: \"", nm, "\"")
    )
  }
})

test_that("get_priors shows one row per element of each parameter", {
  weights <- logistic_normal_pr(
    inv_gamma_pr(3, c(1, 2)),
    tilt = list(R = 3, scale = 2)
  )
  x <- structure(
    list(
      priors = list(
        intercept = normal_pr(c(0, 1), 1),
        sigma = half_normal_pr(c(1, 2)),
        weights = weights
      ),
      setup = list(treated = c("A", "B"), unit = "id")
    ),
    class = "bscmfit"
  )
  tilt <- ", tilt = list(R = 3, scale = 2))"
  expect_equal(
    get_priors(x),
    data.frame(
      parameter = c(
        "intercept", "intercept", "sigma", "sigma", "weights", "weights",
        "concentration", "concentration", "tilt"
      ),
      id = c(rep(c("A", "B"), 4), NA),
      prior = c(
        "normal(0, 1)",
        "normal(1, 1)",
        "half_normal(1)",
        "half_normal(2)",
        paste0("logistic_normal(concentration = inv_gamma(3, 1)", tilt),
        paste0("logistic_normal(concentration = inv_gamma(3, 2)", tilt),
        "inv_gamma(3, 1)",
        "inv_gamma(3, 2)",
        "normal(0, 2)"
      )
    )
  )
  x$priors$weights <- logistic_normal_pr(c(2, 2))
  x$priors$coef <- normal_pr(0, c(1, 2))
  x$setup$beta_names <- c("x", "z")
  d <- get_priors(x)
  expect_named(d, c("parameter", "term", "id", "prior"))
  expect_equal(
    d$parameter,
    c(rep(c("intercept", "sigma", "weights"), each = 2), "coef", "coef")
  )
  expect_equal(d$term, c(rep(NA, 6), "x", "z"))
})
