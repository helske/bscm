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
  omega_prior,
  X_y = NULL,
  X_z = NULL,
  tv_idx = NULL,
  df = 0,
  cv = 0L,
  prior_only = FALSE
) {
  N <- ncol(Y)
  T_total <- nrow(Y)
  J <- ncol(Z)
  standata <- list(
    "T" = T_total,
    T_pre = array(T_pre),
    N = N,
    J = J,
    y = t(Y),
    Z = Z,
    pr_rate_sigma = array(1 / x$sd_e),
    kappa = omega_prior$kappa,
    dirichlet_omega = omega_prior$distribution == "dirichlet",
    cv = cv,
    likelihood = !prior_only
  )
  if (icpt) {
    standata$pr_mean_intercept <- array(x$mean_y)
    standata$pr_sd_intercept <- array(2 * x$sd_e)
  }
  if (!is.null(X_y)) {
    K <- dim(X_y)[3]
    sd_ex <- x$md_sd_e / x$md_sd_x
    standata <- c(
      standata,
      list(
        K = K,
        X_y = X_y,
        X_z = X_z,
        pr_mean_beta = array(0, K),
        pr_sd_beta = array(2 * sd_ex)
      )
    )
    if (!is.null(tv_idx)) {
      L <- length(tv_idx)
      pr_rate_sigma_gamma <- array(2 / sd_ex[tv_idx] * sqrt(df))

      B <- splines::bs(seq_len(T_total), df = df, intercept = TRUE)
      d <- diag(df)
      M <- B %*% lower.tri(d, diag = TRUE)
      T0 <- min(T_pre)
      a <- c(crossprod(M, rep(1:0, c(T0, T_total - T0))))
      Q <- qr.Q(qr(cbind(a, d)))[, -1]
      A <- M %*% Q

      standata <- c(
        standata,
        list(
          L = L,
          tv_idx = array(tv_idx),
          pr_rate_sigma_gamma = array(pr_rate_sigma_gamma),
          D = df,
          A = A
        )
      )
    }
  }
  standata
}

create_inits <- function(x, omega_prior, error = "iid") {
  inits <- list(
    sigma = array(stats::runif(x$N, 0.9, 1.1) / x$pr_rate_sigma)
  )
  if (error == "ar1") {
    inits$rho <- array(stats::runif(x$N, 0, 0.5))
  }
  if (omega_prior$distribution == "logistic_normal") {
    inits$eta <- matrix(0, x$N, x$J - 1)
    inits$omega_ <- matrix(1 / x$J, 0, x$J)
  } else {
    inits$eta <- matrix(0, 0, x$J - 1)
    inits$omega_ <- matrix(1 / x$J, x$N, x$J)
  }
  if (!is.null(x$pr_mean_intercept)) {
    inits$a <- array(
      stats::rnorm(x$N, x$pr_mean_intercept, 0.1 * x$pr_sd_intercept)
    )
  }
  if (!is.null(x$pr_mean_beta)) {
    inits$beta <- array(stats::rnorm(x$K, x$pr_mean_beta, 0.1 * x$pr_sd_beta))
    if (!is.null(x$pr_rate_sigma_gamma)) {
      inits$xi <- matrix(0, x[["D"]] - 1, x$L)
      inits$sigma_gamma <- array(
        stats::runif(x$L, 1, 3) / x$pr_rate_sigma_gamma
      )
    }
  }
  inits
}
