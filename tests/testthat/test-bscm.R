test_that("Incorrect arguments to bscm result in meaningful error message", {
  expect_error(
    bscm(),
    "Argument `formula` is missing\\."
  )
  expect_error(
    bscm(r ~ 1, data = single_treated, treatment = "treatment"),
    "Can't find outcome variable `r` in `data`\\."
  )
  expect_error(
    bscm(~1, data = single_treated, treatment = "treatment"),
    "Argument `formula` must be a <formula> object with an outcome variable"
  )
  expect_error(
    bscm(c(y, x) ~ 1, data = single_treated, treatment = "treatment"),
    "Argument `formula` must be a <formula> object with one outcome variable."
  )
  expect_error(
    bscm(y ~ 1),
    "Argument `data` is missing\\."
  )
  expect_error(
    bscm(y ~ 1, data = 5),
    "Argument `data` must be a <data\\.frame> object\\."
  )
  expect_error(
    bscm(y ~ 1, data = single_treated, treatment = "z"),
    "Can't find treatment variable `z` in `data`\\."
  )
  expect_error(
    bscm(y ~ 1, data = single_treated, treatment = "x"),
    "Variable `x` in `data` should contain either logical or binary values"
  )
  expect_error(
    bscm(y ~ 1, data = single_treated[, -1], treatment = "treatment"),
    "Can't find time index variable `time` in `data`\\."
  )
  expect_error(
    bscm(y ~ 1, data = single_treated, treatment = "treatment", time = "t"),
    "Can't find time index variable `t` in `data`\\."
  )
  expect_error(
    bscm(y ~ 1, data = single_treated[, -2], treatment = "treatment"),
    "Can't find unit index variable `id` in `data`\\."
  )
})


test_that("compute_descriptives() computes correct summary statistics", {
  set.seed(1)
  T_total <- 15L
  T_pre <- c(10L, 8L)
  N <- 2L
  J <- 3L
  K <- 2L
  Y <- matrix(rnorm(T_total * N), T_total, N)
  Z <- matrix(rnorm(T_total * J), T_total, J)
  X <- simplify2array(
    lapply(seq_len(K), \(k) matrix(rnorm((N + J) * T_total), N + J, T_total))
  )
  X_y <- X[seq_len(N), , ]
  X_z <- X[N + seq_len(J), , ]
  x <- compute_descriptives(Y, Z, T_pre, X_y, X_z)

  expect_equal(x$mean_y[1], mean(Y[seq_len(T_pre[1]), 1]))
  expect_equal(x$mean_y[2], mean(Y[seq_len(T_pre[2]), 2]))
  expect_equal(
    x$sd_e[1],
    min(
      sd(Y[seq_len(T_pre[1]), 1]),
      sd(Y[seq_len(T_pre[1]), 1] - rowMeans(Z[seq_len(T_pre[1]), ]))
    )
  )
  expect_length(x$sd_x, K)
})
