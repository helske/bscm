#' Compute descriptive statistics to be used with default priors
#' @noRd
compute_descriptives <- function(
  Y,
  Z,
  T_pre,
  X_y = NULL,
  X_z = NULL,
  beta_names = NULL
) {
  pooled_within_sd <- function(s, n) {
    sqrt(stats::weighted.mean(s^2, w = n - 1))
  }
  N <- ncol(Y)
  sd_y_by_unit <- vapply(
    seq_len(N),
    \(i) sd(Y[seq_len(T_pre[i]), i], na.rm = TRUE),
    numeric(1)
  )
  constant_sd <- which(sd_y_by_unit < 1e-8)
  stopifnot_(
    length(constant_sd) == 0L,
    c(
      "Outcome variable cannot be constant in the pre-treatment period.",
      i = "Found SD < {1e-8} for unit {names(constant_sd)}."
    )
  )
  min_T_pre <- min(T_pre)
  sd_z <- apply(Z[seq_len(min_T_pre), , drop = FALSE], 2, sd)
  constant_sd <- which(sd_z < 1e-8)
  stopifnot_(
    length(constant_sd) == 0L,
    c(
      "Outcome variable cannot be constant in the pre-treatment period.",
      i = "Found SD < {1e-8} for unit {names(constant_sd)}."
    )
  )

  mean_y <- vapply(
    seq_len(N),
    \(i) mean(Y[seq_len(T_pre[i]), i], na.rm = TRUE),
    numeric(1)
  )
  # residual SD with uniform donor weights (ignoring predictors)
  mean_z <- rowMeans(Z)
  sd_yz <- vapply(
    seq_len(N),
    \(i) sd(Y[seq_len(T_pre[i]), i] - mean_z[seq_len(T_pre[i])], na.rm = TRUE),
    numeric(1)
  )
  out <- list(
    mean_y = mean_y,
    sd_e = pmin(sd_y_by_unit, sd_yz)
  )
  if (!is.null(X_y)) {
    out$sd_y <- pooled_within_sd(sd_y_by_unit, T_pre)
    K <- dim(X_y)[3]
    sd_x_by_unit <- matrix(nrow = N, ncol = K)
    for (i in seq_len(N)) {
      sd_x_by_unit[i, ] <- apply(
        X_y[i, seq_len(T_pre[i]), , drop = FALSE],
        3,
        sd
      )
    }
    constant_sd <- which(
      stats::setNames(
        apply(sd_x_by_unit, 2, max) < 1e-8,
        beta_names
      )
    )
    stopifnot_(
      length(constant_sd) == 0,
      c(
        "{cli::qty(length(constant_sd))}
    Model has predictor{?s} which {?is/are} constant in
    the pre-treatment period for all treated units.",
        i = "Found {?a/} constant predictor{?s} {names(constant_sd)}."
      )
    )
    out$sd_x <- apply(
      sd_x_by_unit,
      2,
      pooled_within_sd,
      n = T_pre
    )
  }
  out
}

prior_to_stan <- function(x, type) {
  max_npar <- switch(
    type,
    sigma = 2L,
    omega = 1L,
    intercept = 3L,
    beta = 3L,
    kappa = 2L,
    rho = 2L
  )
  if (type == "omega") {
    pars <- array(0, dim = max(0, x$length))
  } else {
    pars <- matrix(0, max(0, x$length), max_npar)
  }

  if (is.null(x)) {
    dist <- 0
  } else {
    dist <- c(
      normal = 1L,
      student_t = 2L,
      exponential = 3L,
      gamma = 4L,
      beta = 5L,
      half_normal = 6L,
      dirichlet = 7L,
      logistic_normal = 8L
    )[x$distribution]
    switch(
      x$distribution,
      normal = {
        pars[, 1] <- x$location
        pars[, 2] <- x$scale
      },
      student_t = {
        pars[, 1] <- x$df
        pars[, 2] <- x$location
        pars[, 3] <- x$scale
      },
      gamma = {
        pars[, 1] <- x$shape
        pars[, 2] <- x$rate
      },
      beta = {
        pars[, 1] <- x$shape1
        pars[, 2] <- x$shape2
      },
      exponential = {
        pars[, 1] <- x$rate
      },
      half_normal = {
        pars[, 1] <- x$scale
      },
      dirichlet = {
        pars[] <- x$concentration
      },
      logistic_normal = {
        pars[] <- x$scale
      }
    )
  }
  stats::setNames(
    list(dist, pars),
    paste0("pr_", c("dist_", "pars_"), type)
  )
}
create_standata <- function(
  setup,
  priors,
  Y,
  Z,
  X_y = NULL,
  X_z = NULL,
  spline_def = NULL,
  prior_only = FALSE
) {
  N <- ncol(Y)
  T_total <- nrow(Y)
  J <- ncol(Z)
  K <- length(setup$beta_names)
  L <- length(setup$gamma_names)
  if (K == 0) {
    X_y <- array(0, c(0, T_total, K))
    X_z <- array(0, c(0, T_total, K))
  }
  D <- spline_def$D %||% 0
  A <- spline_def$A %||% matrix(0, 0, 0)
  noncentered_xi <- as.integer(spline_def$noncentered %||% 0)

  c(
    list(
      use_alpha = as.integer(setup$has_icpt),
      use_beta = as.integer(setup$has_x),
      use_gamma = as.integer(setup$has_w),
      use_ar1 = as.integer(setup$has_ar1),
      J = J,
      N = N,
      "T" = T_total,
      T_pre = array(setup$T_pre),
      K = K,
      L = L,
      D = D,
      y = t(Y),
      Z = Z,
      X_y = X_y,
      X_z = X_z,
      tv_idx = array(setup$tv_idx),
      A = A,
      likelihood = as.integer(!prior_only),
      noncentered_xi = noncentered_xi,
      missing_idx = setup$missing_idx,
      n_missing = nrow(setup$missing_idx)
    ),
    unlist(
      lapply(
        c("sigma", "omega", "intercept", "beta", "kappa", "rho"),
        \(i) prior_to_stan(priors[[i]], i)
      ),
      recursive = FALSE
    )
  )
}
# x: standata, d: descriptives
create_inits <- function(x, d, spline_def) {
  eta <- matrix(0, x$N, x$J - 1)
  sigma <- array(stats::runif(x$N, 0.5, 1.5) * d$sd_e)
  a <- beta <- rho <- kappa <- array(0, 0)
  xi <- matrix(0, 0, 0)
  if (x$use_alpha) {
    a <- array(stats::rnorm(x$N, d$mean_y, 0.25 * d$sd_e))
  }
  if (x$use_beta) {
    sd_yx <- d$sd_y / d$sd_x
    beta <- array(stats::rnorm(x$K, 0, 0.1 * sd_yx))
    if (x$use_gamma) {
      xi <- matrix(0, x$D, x$L)
      kappa <- array(stats::runif(x$L, 0.01, 0.05))
    }
  }
  if (x$use_ar1) {
    rho <- array(stats::runif(x$N, 0, 0.5))
  }
  dplyr::lst(eta, sigma, a, beta, xi, kappa, rho)
}
