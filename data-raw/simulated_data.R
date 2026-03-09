## code to used to create `simulated_data` dataset
library(dplyr)
set.seed(454546)
# number of donor series
J <- 30
# total number of time periods, last 10 correspond for post-treatment
T_total <- 50
# two factors
psi <- cbind(
  c(1 + arima.sim(list(ar = 0.3), n = T_total)),
  c(2 + arima.sim(list(ar = 0.9), n = T_total, sd = 0.5))
)
lambda <- rbind(
  rnorm(J + 1, 2, 0.5),
  runif(J + 1)
)
y <- c(psi %*% lambda)

# additional predictor
x <- 2 * rep(colSums(lambda), each = T_total) + 
  psi[, 1, drop = FALSE] %*% lambda[1, , drop = FALSE] + 
  rnorm(length(y), 0, 2)
simulated_data <- data.frame(
  time = seq_len(T_total) - (T_total - 9),
  id = rep(seq_len(J + 1), each = T_total),
  y = y,
  x = c(x),
  beta = 1,
  psi1 = psi[, 1],
  psi2 = psi[, 2],
  lambda1 = rep(lambda[1, ], each = T_total),
  lambda2 = rep(lambda[2, ], each = T_total)
) |> 
  dplyr::mutate(
    treatment = ifelse(id == 1 & time >= 0, 1, 0),
    .before = beta
  ) |> 
  dplyr::mutate(
    # true effect is log(1 + t), t = time points since treatment
    tau = ifelse(treatment, log(1 + cumsum(treatment)), 0),
    .by = id,
    .before = beta
  ) |> 
  dplyr::mutate(
    y = y + beta * x + tau * treatment + rnorm(n())
  )
usethis::use_data(simulated_data, overwrite = TRUE)
