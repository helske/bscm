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
