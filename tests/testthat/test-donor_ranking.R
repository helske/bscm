test_that("donor_ranking returns a character vector", {
  d <- donor_ranking(fit1_int)
  expect_type(d, "character")
})

test_that("donor_ranking length equals number of donors", {
  d <- donor_ranking(fit1_int)
  expect_length(d, get_J(fit1_int))
})

test_that("donor_ranking contains exactly the donor IDs", {
  d <- donor_ranking(fit1_int)
  expect_setequal(d, get_donors(fit1_int))
})

test_that("donor_ranking returns donors in descending weight order", {
  d <- donor_ranking(fit1_int)
  w <- donor_weights(fit1_int, probs = numeric(0))
  w_ordered <- w |> dplyr::arrange(dplyr::desc(mean))
  expect_equal(d, w_ordered$donor)
})

test_that("donor_ranking errors for non-bscmfit non-vsel input", {
  expect_error(
    donor_ranking(list()),
    "must be an object of class"
  )
  expect_error(
    donor_ranking("fit1_int"),
    "must be an object of class"
  )
})
