## code to used to create `simulated_example` dataset
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
  c(2 + arima.sim(list(ar = 0.95, sd = 0.5), n = T_total)),
  log(4 + seq_len(T_total)) + 0.2 * abs(rnorm(T_total))
)
lambda <- t(MCMCpack::rdirichlet(J, rep(0.8, 3)))
w <- c(MCMCpack::rdirichlet(1, rep(0.15, J)))
# enforce convex hull for treated unit _before_ adding noise
lambda <- cbind(lambda %*% w, lambda)
y <- psi %*% lambda

# common time-varying intercept
alpha <- cumsum(cumsum(rnorm(T_total, 0, 0.05)))
# additional predictors:
# time-invariant predictor correlated with the loadings of the first factor
x1 <- rep(rnorm(J + 1, lambda[1, ], 0.5), each = T_total)
# time-varying predictor correlated with second factor
x2 <- 0.3 * psi[, 2] + rnorm(T_total * (J + 1))
# time-varying effect for x1
beta1 <- 0.5 * psi[, 3] + rnorm(T_total)
beta2 <- 0.5 # constant effect for x2

simulated_example <- data.frame(
  time = seq_len(T_total) - (T_total - 9),
  id = paste0("id_", rep(seq_len(J + 1), each = T_total)),
  x1 = x1,
  x2 = x2,
  y = c(y),
  alpha = alpha,
  beta1 = beta1,
  beta2 = beta2
) |> 
  dplyr::mutate(
    id = factor(id, levels = paste0("id_", seq_len(J + 1))),
    treatment = ifelse(id == "id_1" & time >= 0, 1, 0)
  ) |> 
  dplyr::mutate(
    # true effect is log(1 + t), t = time points since treatment
    tau = ifelse(treatment, log(1 + cumsum(treatment)), 0),
    .by = id
  ) |> 
  dplyr::mutate(
    y = y + alpha + beta1 * x1 + beta2 * x2 + tau + rnorm(n(), 0, 0.25)
  )
usethis::use_data(simulated_example, overwrite = TRUE)
