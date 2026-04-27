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
    bscm( ~ 1, data = single_treated, treatment = "treatment"),
    "Argument `formula` must be a <formula> object with an outcome variable on the left-hand side\\."
  )
  expect_error(
    bscm(c(y, x) ~ 1, data = single_treated, treatment = "treatment"),
    "Argument `formula` must be a <formula> object with one outcome variable on the left-hand side\\."
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
    bscm(y ~ 1, data = single_treated),
    "Argument `treatment` is missing\\."
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


test_that("bscm_stats() computes correct summary statistics", {
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
  x <- bscm_stats(Y, Z, T_pre, X)
  
  expect_equal(x$mean_y[1], mean(Y[seq_len(T_pre[1]), 1]))
  expect_equal(x$mean_y[2], mean(Y[seq_len(T_pre[2]), 2]))
  expect_equal(
    x$sd_e[1], 
    sd(Y[seq_len(T_pre[1]), 1] - rowMeans(Z[seq_len(T_pre[1]), ]))
  )
  expect_true(all(x$sd_e >= 1))
  expect_gte(x$md_sd_e, 1)
  expect_length(x$md_sd_x, K)
  expect_true(all(x$md_sd_x >= 1))
})

test_that("create_standata() returns correct list structure", {
  set.seed(1)
  kappa <- 2
  T_total <- 15L
  T_pre <- c(10L, 8L)
  N <- 2L
  J <- 3L
  Y <- matrix(rnorm(T_total * N), T_total, N)
  Z <- matrix(rnorm(T_total * J), T_total, J)
  x <- bscm_stats(Y, Z, T_pre)
  d <- create_standata(x, T_pre, Y, Z, icpt = TRUE, kappa)
 
  expect_named(
    d,
    c("T", "T_pre", "N", "J", "y", "Z",
      "pr_rate_sigma", "kappa",  "pr_mean_intercept", "pr_sd_intercept")
  )
  expect_equal(d$T, T_total)
  expect_equal(d$N, N)
  expect_equal(d$J, J)
  expect_equal(d$kappa, kappa)
  expect_equal(d$pr_mean_intercept, array(x$mean_y))
  expect_equal(d$pr_rate_sigma, array(1 / x$sd_e))
})

test_that("create_inits() returns correct list structure", {
  set.seed(1)
  kappa <- 2
  T_total <- 15L
  T_pre <- c(10L, 8L)
  N <- 2L
  J <- 3L
  Y <- matrix(rnorm(T_total * N), T_total, N)
  Z <- matrix(rnorm(T_total * J), T_total, J)
  x <- bscm_stats(Y, Z, T_pre)
  d <- create_standata(x, T_pre, Y, Z, icpt = FALSE, kappa)
  inits <- create_inits(d, logistic_normal())
  expect_named(inits, c("eta", "sigma"))
  expect_equal(dim(inits$eta), c(N, J))
  expect_length(inits$sigma, N)
  d <- create_standata(x, T_pre, Y, Z, icpt = TRUE, kappa)
  inits <- create_inits(d, logistic_normal())
  expect_named(inits, c("eta", "sigma", "a"))
  expect_equal(dim(inits$eta), c(N, J))
  expect_length(inits$sigma, N)
  expect_length(inits$a, N)
  inits_dir <- create_inits(d, dirichlet(kappa = 1))
  expect_named(inits_dir, c("omega", "sigma", "a"))
  expect_length(inits_dir$omega, N)
  expect_equal(inits_dir$omega[[1]], rep(1 / J, J))
})
