bscm_stats <- function(Y, Z, T_pre, X = NULL) {
  N <- ncol(Y)
  sd_y <- vapply(seq_len(N), \(i) sd(Y[seq_len(T_pre[i]), i]), numeric(1))
  stopifnot_(
    all(sd_y > sqrt(.Machine$double.eps)),
    "Outcome variable cannot be constant in the pre-treatment period. 
    Found `sd(y) < sqrt(.Machine$double.eps)`."
  )
  min_T_pre <- min(T_pre)
  sd_z <- apply(Z[seq_len(min_T_pre), , drop = FALSE], 2, sd)
  stopifnot_(
    any(sd_z > sqrt(.Machine$double.eps)),
    "Outcome variable cannot be constant in the pre-treatment period. 
    Found `sd(z) < sqrt(.Machine$double.eps)`."
  )

  mean_y <- vapply(seq_len(N), \(i) mean(Y[seq_len(T_pre[i]), i]), numeric(1))
  # residual SD with uniform donor weights (ignoring predictors)
  mean_sc <- rowMeans(Z)
  sd_e <- vapply(
    seq_len(N),
    \(i) {
      sd(Y[seq_len(T_pre[i]), i] - mean_sc[seq_len(T_pre[i])])
    },
    numeric(1)
  )
  out <- list(
    mean_y = mean_y,
    sd_e = pmax(1, sd_e),
    md_sd_e = pmax(1, stats::median(sd_e))
  )
  if (!is.null(X)) {
    sd_x_by_unit <- apply(X[, seq_len(min_T_pre), , drop = FALSE], c(1, 3), sd)
    out$md_sd_x <- pmax(1, apply(sd_x_by_unit, 2, stats::median))
  }
  out
}

create_standata <- function(
  x,
  T_pre,
  Y,
  Z,
  icpt,
  kappa,
  X_y = NULL,
  X_z = NULL,
  tv_idx = NULL
) {
  N <- ncol(Y)
  T_total <- nrow(Y)
  J <- ncol(Z)
  standata <- list(
    T = T_total,
    T_pre = array(T_pre),
    N = N,
    J = J,
    y = t(Y),
    Z = Z,
    pr_rate_sigma = array(1 / x$sd_e),
    kappa = kappa
  )
  if (icpt) {
    standata$pr_mean_intercept <- array(x$mean_y)
    standata$pr_sd_intercept <- array(2 * x$sd_e)
  }
  if (!is.null(X_y)) {
    K <- dim(X_y)[3]
    pr_sd_beta <- 2 * x$md_sd_e / x$md_sd_x
    standata <- c(
      standata,
      list(
        K = K,
        X_y = X_y,
        X_z = aperm(X_z, c(3, 2, 1)),
        pr_mean_beta = array(0, K),
        pr_sd_beta = array(pr_sd_beta)
      )
    )
    if (!is.null(tv_idx)) {
      L <- length(tv_idx)
      # mode at 0.1sigma, P(sigma_gamma < 0.5sigma) = 0.96
      pr_rate_sigma_gamma <- rep(10 / x$md_sd_e, L)
      standata <- c(
        standata,
        list(
          L = L,
          tv_idx = array(tv_idx),
          pr_rate_sigma_gamma = array(pr_rate_sigma_gamma)
        )
      )
    }
  }
  standata
}

create_inits <- function(x, omega_prior) {
  inits <- list(
    sigma = array(stats::runif(x$N, 0.9, 1.1) / x$pr_rate_sigma)
  )
  if (omega_prior$distribution == "logistic_normal") {
    inits$eta <- matrix(0, x$N, x$J - 1)
  } else {
    inits$omega <- matrix(1 / x$J, x$N, x$J)
  }
  if (!is.null(x$pr_mean_intercept)) {
    inits$a <- array(
      stats::rnorm(x$N, x$pr_mean_intercept, 0.5 * x$pr_sd_intercept)
    )
  }
  if (!is.null(x$pr_mean_beta)) {
    inits$beta <- array(stats::rnorm(x$K, x$pr_mean_beta, 0.1 * x$pr_sd_beta))
    if (!is.null(x$pr_rate_sigma_gamma)) {
      inits$gamma_raw <- matrix(0, x$T, x$L)
      inits$sigma_gamma <- array(
        stats::runif(x$L, 0.9, 1.1) / x$pr_rate_sigma_gamma
      )
    }
  }
  inits
}
