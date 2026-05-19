## code to used to create `single_treated` dataset
library(dplyr)
set.seed(454546)
# number of donor series
J <- 30
# total number of time periods, last 10 correspond for post-treatment
T_total <- 40
# two factors
psi <- cbind(
  c(-1 + arima.sim(list(ar = 0.9), n = T_total, sd = sqrt(0.19))),
  0.5 * (-sqrt(seq_len(T_total)) + runif(T_total))
)
lambda <- rbind(
  rgamma(J, 2, 2),
  rbeta(J, 4, 2)
)
lambda <- cbind(lambda %*% rep(1 / J, J), lambda)
alpha <- rep(rnorm(J + 1, 0, 2), each = T_total)
y <- alpha + psi %*% lambda

# additional predictor
x <- matrix(0, T_total, J + 1)
for (j in seq_len(J + 1)) {
  x[, j] <- sin(0.5 * seq_len(T_total)) + 0.2 * cumsum(rnorm(T_total))
}
single_treated <- data.frame(
  time = seq_len(T_total) - (T_total - 9),
  id = rep(seq_len(J + 1), each = T_total),
  y = c(y),
  x = c(x),
  alpha = alpha,
  psi1 = psi[, 1],
  psi2 = psi[, 2],
  lambda1 = rep(lambda[1, ], each = T_total),
  lambda2 = rep(lambda[2, ], each = T_total)
) |>
  dplyr::mutate(
    treatment = ifelse(id == 1 & time >= 0, 1, 0),
    .before = alpha
  ) |>
  dplyr::mutate(
    y = 8 + y + x + (time + 1) * treatment + rnorm(n(), 0, 0.5)
  )
usethis::use_data(single_treated, overwrite = TRUE)
