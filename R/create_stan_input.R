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
    sd_e = sd_e,
    md_sd_e = stats::median(sd_e)
  )
  if (!is.null(X)) {
    sd_x_by_unit <- apply(X[, seq_len(min_T_pre), , drop = FALSE], c(1, 3), sd)
    out$md_sd_x <- apply(sd_x_by_unit, 2, stats::median)
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
  tv_idx = numeric(0),
  df = 0,
  ar1_error = FALSE,
  prior_only = FALSE
) {
  N <- ncol(Y)
  T_total <- nrow(Y)
  J <- ncol(Z)
  K <- 0
  L <- length(tv_idx)
  D <- 1
  if (L == 0) {
    A <- matrix(0, T_total, 0)
    pr_sd_sigma_gamma <- array(0)
  }
  if (is.null(X_y)) {
    X_y <- array(0, c(N, T_total, 0))
    X_z <- array(0, c(J, T_total, 0))
    sd_ex <- 0
  } else {
    K <- dim(X_y)[3]
    sd_ex <- x$md_sd_e / x$md_sd_x
    if (L > 0) {
      D <- df
      pr_sd_sigma_gamma <- array(sd_ex[tv_idx] * sqrt(T_total / D))
      B <- splines::bs(seq_len(T_total), df = df, intercept = TRUE)
      d <- diag(df)
      M <- B %*% lower.tri(d, diag = TRUE)
      T0 <- min(T_pre)
      a <- c(crossprod(M, rep(1:0, c(T0, T_total - T0))))
      Q <- qr.Q(qr(cbind(a, d)))[, -1]
      A <- M %*% Q
    }
  }

  list(
    J = J,
    N = N,
    "T" = T_total,
    T_pre = array(T_pre),
    K = K,
    L = L,
    D = D,
    icpt = as.integer(icpt),
    ar1 = as.integer(ar1_error),
    dirichlet_omega = as.integer(omega_prior$distribution == "dirichlet"),
    kappa = omega_prior$kappa,
    y = t(Y),
    Z = Z,
    X_y = X_y,
    X_z = X_z,
    tv_idx = array(tv_idx),
    A = A,
    pr_rate_sigma = array(1 / x$sd_e),
    pr_mean_intercept = array(x$mean_y, N * icpt),
    pr_sd_intercept = array(2 * x$sd_e, N * icpt),
    pr_mean_beta = array(0, K),
    pr_sd_beta = array(2 * sd_ex, K),
    pr_sd_sigma_gamma = pr_sd_sigma_gamma,
    pr_shape1_rho = array(2, N * ar1_error),
    pr_shape2_rho = array(2, N * ar1_error),
    likelihood = as.integer(!prior_only)
  )
}

create_inits <- function(x) {
  if (x$dirichlet_omega) {
    omega_ <- matrix(1 / x$J, x$N, x$J)
    eta <- matrix(0, 0, x$J - 1)
  } else {
    omega_ <- matrix(1 / x$J, 0, x$J)
    eta <- matrix(0, x$N, x$J - 1)
  }
  sigma <- array(stats::runif(x$N, 0.9, 1.1) / x$pr_rate_sigma)
  a <- beta <- rho <- sigma_gamma <- array(0, 0)
  xi <- matrix(0, 0, 0)
  if (x$icpt) {
    a <- array(stats::rnorm(x$N, x$pr_mean_intercept, 0.1 * x$pr_sd_intercept))
  }
  if (x$K > 0) {
    beta <- array(stats::rnorm(x$K, x$pr_mean_beta, 0.1 * x$pr_sd_beta))
    if (x$L > 0) {
      xi <- matrix(0, x$D - 1, x$L)
      #sigma_gamma <- array(stats::runif(x$L, 1, 3) / x$pr_rate_sigma_gamma)
      sigma_gamma <- array(stats::runif(x$L, 0.1, 0.5) * x$pr_sd_sigma_gamma)
    }
  }
  if (x$ar1) {
    rho <- array(stats::runif(x$N, 0, 0.5))
  }
  dplyr::lst(omega_, eta, sigma, a, beta, xi, sigma_gamma, rho)
}
