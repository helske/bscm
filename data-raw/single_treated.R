## code to used to create `single_treated` dataset

set.seed(454546)
T_0 <- 30
T_post <- 5
T_total <- T_0 + T_post
J <- 50
psi <- cbind(
  4 + rnorm(T_total),
  2 * sin(0.1 * seq_len(T_total)) + 0.2 * rnorm(T_total),
  0.1 * seq_len(T_total) + rnorm(T_total)
)
# factor loadings
lambda <- rbind(
  rnorm(J),
  rnorm(J),
  2 + rnorm(J)
)
lambda <- cbind(c(0, 0.5, 1), lambda)
f <- psi %*% lambda
# unit- and time-specific intercepts
alpha <- rnorm(J + 1)
delta <- rnorm(T_total)

# N(0, 0.5) noise term
epsilon <- rnorm((J + 1) * T_total, 0, 0.5)

# time-varying covariate
x <- matrix(0, T_total, J + 1)
for (j in seq_len(J + 1)) {
  x[, j] <- lambda[3, j] * psi[, 3] + arima.sim(list(ar = 0.9), T_total)
}

single_treated <- data.frame(
  time = seq_len(T_total) - T_0 - 1,
  id = rep(seq_len(J + 1), each = T_total),
  x = c(x),
  alpha = rep(alpha, each = T_total), 
  delta = delta,
  f = c(f),
  epsilon = epsilon
) |> 
  dplyr::mutate(
    treatment = as.integer(id == 1 & time >= 0),
    y = alpha + delta + f + (1 + time) * treatment + epsilon,
  ) |> 
  dplyr::select(time, id, y, x, treatment)

usethis::use_data(single_treated, overwrite = TRUE)
