compute_ess_dirichlet <- function(kappa, J, nsim) {
  x <- matrix(rgamma(nsim * J, kappa, 1), J, nsim)
  colSums(x)^2 / colSums(x^2)
}
compute_ess_logistic_normal <- function(x) {
  w <- apply(x, 2, \(x) {
    expx <- exp(x - max(x))
    expx / sum(expx)
  })
  1 / colSums(w^2)
}

# number of donors
J_seq <- seq(10, 500, by = 10)
# number of Monte Carlo replications
nsim <- 1e4

# scale parameter of logistic normal
kappa_seq <- seq(0.1, 5, by = 0.1)
K <- length(kappa_seq)

set.seed(23546)
results_logistic_normal <-
  lapply(
    J_seq,
    \(J) {
      qs <- matrix(NA, 3, K)
      x <- matrix(rnorm(nsim * J), J, nsim)
      for (i in seq_along(kappa_seq)) {
        ess <- compute_ess_logistic_normal(kappa_seq[i] * x)
        qs[, i] <- quantile(ess, c(0.05, 0.5, 0.95))
      }
      data.frame(
        J = J,
        kappa = kappa_seq,
        median_ess = qs[2, ],
        q5_ess = qs[1, ],
        q95_ess = qs[3, ]
      )
    }
  ) |>
  bind_rows()

# values of concentration parameter of symmetric Dirichlet distribution
kappa_seq <- c(seq(0.1, 1, by = 0.05), seq(1.5, 10, by = 0.5))
K <- length(kappa_seq)

set.seed(23546)
results_dirichlet <- lapply(
  J_seq,
  \(J) {
    qs <- matrix(NA, 3, K)
    for (i in seq_along(kappa_seq)) {
      ess <- compute_ess_dirichlet(kappa_seq[i], J, nsim)
      qs[, i] <- quantile(ess, c(0.05, 0.5, 0.95))
    }
    data.frame(
      J = J,
      kappa = kappa_seq,
      median_ess = qs[2, ],
      q5_ess = qs[1, ],
      q95_ess = qs[3, ]
    )
  }
) |>
  bind_rows()

kappa_lookup <- bind_rows(
  dirichlet = results_dirichlet,
  logistic_normal = results_logistic_normal,
  .id = "distribution"
) |>
  mutate(
    median_r_ess = median_ess / J,
    q5_r_ess = q5_ess / J,
    q95_r_ess = q95_ess / J
  )
usethis::use_data(kappa_lookup)
