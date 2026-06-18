## code to used to create `multiple_treated` dataset
set.seed(235)
# number of donor series
J <- 50
# number of treated series
N <- 3
# total number of time periods, last 5 correspond for post-treatment
T_total <- 40
# two factors
psi <- cbind(
  c(arima.sim(list(ar = 0.9), n = T_total)),
  c(-0.2 * seq_len(T_total) + arima.sim(list(ar = 0.3), n = T_total))
)
lambda <- rbind(
  rnorm(J + N, 1, 1),
  rbeta(J + N, 2, 2)
)
lambda[2, 1:N] <- rbeta(N, 5, 2)
f <- c(psi %*% lambda)
# additional predictors
x <- matrix(0, T_total, J + N)
for (j in seq_len(J + N)) {
  x[, j] <- sin(0.5 * seq_len(T_total)) + 0.2 * cumsum(rnorm(T_total))
}
z <- matrix(0, T_total, J + N)
for (j in seq_len(J + N)) {
  z[, j] <- lambda[1, j] * psi[, 1] + 2 * lambda[2, j] + rnorm(T_total)
}
# time-varying regression coefficient for z
gamma <- cumsum(cumsum(rnorm(T_total, 0, 0.05)))

# noise term
epsilon <- rnorm((J + N) * T_total)

multiple_treated <- data.frame(
  time = seq_len(T_total) - (T_total - 4),
  id = rep(seq_len(J + N), each = T_total),
  f = c(f),
  x = c(x),
  z = c(z),
  gamma = gamma,
  psi1 = psi[, 1],
  psi2 = psi[, 2],
  lambda1 = rep(lambda[1, ], each = T_total),
  lambda2 = rep(lambda[2, ], each = T_total),
  epsilon = epsilon
) |>
  dplyr::mutate(
    treatment = as.integer(id <= N & time >= 0),
    y = 5 + f + x + gamma * z + 2 * treatment + epsilon
  ) |> 
  dplyr::select(time, id, y, x, z, treatment)
usethis::use_data(multiple_treated, overwrite = TRUE)
