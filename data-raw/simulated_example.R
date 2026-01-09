## code to prepare `simulated_example` dataset goes here
library(MCMCpack)
library(dplyr)
set.seed(454546)
# number of donor series
J <- 30
# total number of time periods, last 10 correspond for post-treatment
T_total <- 50
# three latent factors: N(5, 2), AR(1), and a noisy nonlinear trend
psi <- cbind(
  rnorm(T_total, 5, 1),
  2 + arima.sim(list(ar = 0.9), n = T_total),
  log(10 + seq_len(T_total)) + rnorm(T_total)
)
lambda <- t(MCMCpack::rdirichlet(J, rep(0.8, 3)))
w <- c(MCMCpack::rdirichlet(1, rep(0.15, J)))
# enforce convex hull for treated unit _before_ adding noise
lambda <- cbind(lambda %*% w, lambda)
y <- psi %*% lambda

simulated_example <- data.frame(
  time = seq_len(T_total) - (T_total - 9),
  id = rep(seq_len(J + 1), each = T_total),
  y = c(y)
) |> 
  dplyr::mutate(
    treatment = ifelse(id == 1 & time >= 0, 1, 0)
  ) |> 
  dplyr::mutate(
    # true effect is log(1 + t), t = time points since treatment
    tau = ifelse(treatment, log(1 + cumsum(treatment)), 0),
    .by = id
  ) |> 
  dplyr::mutate(
    y = y + rnorm(n(), 0, 0.25) # add noise
  )
usethis::use_data(simulated_example, overwrite = TRUE)
