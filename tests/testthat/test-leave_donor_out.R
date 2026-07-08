test_that("leave_donor_out validates cumulative argument", {
  expect_error(
    leave_donor_out(fit1_int, cumulative = "yes"),
    "Argument `cumulative` must be either TRUE or FALSE"
  )
})

test_that("leave_donor_out validates probs argument", {
  expect_error(
    leave_donor_out(fit1_int, probs = "a"),
    "Argument `probs` must be a <numeric> vector"
  )
})

test_that("leave_donor_out errors when no saved data", {
  fit_nosave <- fit1_int
  fit_nosave$data <- NULL
  expect_error(
    leave_donor_out(fit_nosave),
    "does not contain the original data"
  )
})
