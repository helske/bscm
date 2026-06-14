test_that("parse_bscm_formula returns formula unchanged when no tv()", {
  f <- y ~ x + z
  parsed <- parse_bscm_formula(f)
  expect_null(parsed$w_formula)
  expect_identical(parsed$w_terms, character(0))
  expect_equal(
    sort(attr(terms(parsed$x_formula), "term.labels")),
    sort(attr(terms(f), "term.labels"))
  )
})

test_that("parse_bscm_formula extracts tv() terms", {
  f <- y ~ x + tv(~z, 7)
  parsed <- parse_bscm_formula(f)
  expect_identical(parsed$w_terms, "z")
  full_terms <- attr(terms(parsed$x_formula), "term.labels")
  expect_true(all(c("x", "z") %in% full_terms))
  w_terms <- attr(terms(parsed$w_formula), "term.labels")
  expect_identical(w_terms, "z")
  expect_equal(attr(terms(parsed$w_formula), "intercept"), 1L)
})

test_that("parse_bscm_formula handles interactions inside tv()", {
  f <- y ~ x * z + tv(~ x + w * r, 10)
  parsed <- parse_bscm_formula(f)
  full_terms <- attr(terms(parsed$x_formula), "term.labels")
  expect_true(all(c("x", "z", "x:z", "w", "r", "w:r") %in% full_terms))
  w_terms <- attr(terms(parsed$w_formula), "term.labels")
  expect_true(all(c("x", "w", "r", "w:r") %in% w_terms))
})

test_that("parse_bscm_formula preserves intercept status", {
  f_int <- y ~ x + tv(~z, 10)
  f_noint <- y ~ 0 + x + tv(~z, 100)
  parsed_int <- parse_bscm_formula(f_int)
  parsed_noint <- parse_bscm_formula(f_noint)
  expect_equal(attr(terms(parsed_int$x_formula), "intercept"), 1L)
  expect_equal(attr(terms(parsed_noint$x_formula), "intercept"), 0L)
  expect_equal(attr(terms(parsed_int$w_formula), "intercept"), 1L)
  expect_equal(attr(terms(parsed_noint$w_formula), "intercept"), 0L)
})

test_that("parse_bscm_formula tv terms are included in full formula", {
  f <- y ~ tv(~ x + z + r * log(w), 6)
  parsed <- parse_bscm_formula(f)
  full_terms <- attr(terms(parsed$x_formula), "term.labels")
  expect_true(all(c("x", "z", "r", "log(w)", "r:log(w)") %in% full_terms))
})

test_that("parse_bscm_formula errors on multiple tv() calls", {
  f <- y ~ tv(~x, 10) + tv(~z)
  expect_error(
    parse_bscm_formula(f),
    "Multiple `tv\\(\\)` terms are not supported"
  )
})

test_that("parse_bscm_formula errors on non-formula tv() argument", {
  f <- y ~ tv("aa")
  expect_error(
    parse_bscm_formula(f),
    "must be a one-sided formula"
  )
})

test_that("parse_bscm_formula errors on two-sided tv() formula", {
  f <- y ~ tv(a ~ x, 5)
  expect_error(
    parse_bscm_formula(f),
    "must be a one-sided formula"
  )
})

test_that("parse_bscm_formula errors on empty tv()", {
  f <- y ~ tv(~ -1)
  expect_error(
    parse_bscm_formula(f),
    "The formula inside `tv\\(\\)` must contain at least one term\\."
  )
  f <- y ~ tv(~1, 9)
  expect_error(
    parse_bscm_formula(f),
    "The formula inside `tv\\(\\)` must contain at least one term"
  )
})
test_that("parse_bscm_formula errors on df < 5", {
  f <- y ~ tv(~x, 2)
  expect_error(
    parse_bscm_formula(y ~ tv(~x, 2)),
    "must be an integer larger than 4"
  )
})

test_that("parse_bscm_formula errors on nested tv()", {
  f <- y ~ tv(~ tv(~x), 9)
  expect_error(
    parse_bscm_formula(f),
    "Nested `tv\\(\\)` terms are not supported"
  )
})
