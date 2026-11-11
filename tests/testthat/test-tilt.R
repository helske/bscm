test_that("create_tilt_data computes standardized directions from R", {
  set.seed(1)
  T_total <- 12L
  J <- 6L
  donors <- letters[seq_len(J)]
  Z <- matrix(rnorm(T_total * J), T_total, J, dimnames = list(NULL, donors))
  x <- create_tilt_data(list(R = 2L, scale = 1.5), Z, c(9L, 8L), donors)
  expect_named(x, c("R", "V", "sd_phi"))
  expect_equal(x$R, 2L)
  expect_equal(x$sd_phi, 1.5)
  expect_equal(dim(x$V), c(J, 2L))
  expect_equal(unname(colMeans(x$V)), c(0, 0))
  expect_equal(unname(apply(x$V, 2, sd)), c(1, 1))
  # directions are based on the first min(T_pre) time points only
  Z2 <- Z
  Z2[9:T_total, ] <- rnorm((T_total - 8L) * J)
  expect_equal(
    create_tilt_data(list(R = 2L, scale = 1.5), Z2, c(9L, 8L), donors),
    x
  )
  # directions are invariant to donor levels and common time effects
  Z3 <- sweep(Z, 2, rnorm(J), "+") + rnorm(T_total)
  x3 <- create_tilt_data(list(R = 2L, scale = 1.5), Z3, c(9L, 8L), donors)
  expect_equal(abs(x3$V), abs(x$V))
})

test_that("create_tilt_data checks the number of directions", {
  Z <- matrix(rnorm(40), 10, 4)
  donors <- letters[1:4]
  expect_equal(
    create_tilt_data(list(R = 3L, scale = 1), Z, 10L, donors)$R,
    3L
  )
  expect_error(
    create_tilt_data(list(R = 4L, scale = 1), Z, 10L, donors),
    "Element \"R\" of `tilt` must be at most 3\\."
  )
  expect_error(
    create_tilt_data(list(R = 3L, scale = 1), Z, 3L, donors),
    "Element \"R\" of `tilt` must be at most 2\\."
  )
})

test_that("create_tilt_data matches the rows of user-defined V to donors", {
  Z <- matrix(rnorm(40), 10, 4)
  donors <- c("a", "b", "c", "d")
  V <- matrix(
    1:10,
    5,
    dimnames = list(c("d", "x", "b", "a", "c"), NULL)
  )
  x <- create_tilt_data(list(R = 2L, V = V, scale = 2), Z, 10L, donors)
  expect_equal(x$R, 2L)
  expect_equal(x$V, unname(V[donors, ]))
  expect_equal(x$sd_phi, 2)
  expect_error(
    create_tilt_data(
      list(R = 2L, V = V[-1, ], scale = 2),
      Z,
      10L,
      donors
    ),
    "Element \"V\" of `tilt` has no rows for donor \"d\"\\."
  )
})

test_that("create_inits adds tilt coefficients only for the tilt model", {
  x <- list(
    N = 2L,
    J = 4L,
    use_alpha = 0L,
    use_beta = 0L,
    use_ar1 = 0L,
    use_zeta = 1L
  )
  d <- list(sd_e = c(1, 2))
  inits <- create_inits(x, d, NULL, marginalize_sigma = TRUE)
  expect_null(inits$phi)
  expect_null(inits$sigma)
  x$R <- 3L
  x$V <- matrix(0, 4, 3)
  inits <- create_inits(x, d, NULL, marginalize_sigma = TRUE)
  expect_equal(dim(inits$phi), c(2L, 3L))
})
