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
    max_T_pre <- max(T_pre)
    out$sd_t <- sqrt(max_T_pre * (max_T_pre + 1) / 12)
  }
  out
}

prior_to_stan <- function(x, type) {

  max_npar <- switch(
    type,
    sigma = 3L,
    omega = 1L,
    alpha = 3L,
    beta = 3L,
    kappa = 3L,
    nu = 3L,
    rho = 2L,
    zeta = 3L
  )
  if (type == "omega") {
    if (x$zeta) x$length <- 0L
    pars <- array(0, dim = max(0, x$length))
  } else {
    pars <- matrix(0, max(0, x$length), max_npar)
  }

  if (is.null(x)) {
    dist <- 0
  } else {
    # codes must match prior_lpdf() and build_omega_lp() in priors.stan,
    # negative codes are used for the simplex priors of omega
    dist <- c(
      beta = 1L,
      cauchy = 2L,
      dirichlet = -1L,
      exponential = 3L,
      gamma = 4L,
      half_normal = 5L,
      inv_gamma = 6L,
      laplace = 7L,
      logistic_laplace = -2L,
      logistic_normal = -3L,
      normal = 8L,
      scaled_inv_chi = 9L,
      sqrt_inv_gamma = 10L,
      student_t = 11L
    )[x$distribution]
    switch(
      x$distribution,
      beta = {
        pars[, 1] <- x$shape1
        pars[, 2] <- x$shape2
      },
      cauchy = {
        pars[, 1] <- x$location
        pars[, 2] <- x$scale
      },
      dirichlet = {
        pars[] <- x$concentration
      },
      exponential = {
        pars[, 1] <- x$rate
      },
      gamma = {
        pars[, 1] <- x$shape
        pars[, 2] <- x$rate
      },
      half_normal = {
        pars[, 1] <- x$scale
      },
      inv_gamma = {
        pars[, 1] <- x$shape
        pars[, 2] <- x$scale
      },
      laplace = {
        pars[, 1] <- x$location
        pars[, 2] <- x$scale
      },
      logistic_laplace = {
        pars[] <- x$concentration
      },
      logistic_normal = {
        pars[] <- x$concentration
      },
      normal = {
        pars[, 1] <- x$location
        pars[, 2] <- x$scale
      },
      scaled_inv_chi = {
        pars[, 1] <- x$df
        pars[, 2] <- x$scale
      },
      sqrt_inv_gamma = {
        pars[, 1] <- x$shape
        pars[, 2] <- x$scale
      },
      student_t = {
        pars[, 1] <- x$df
        pars[, 2] <- x$location
        pars[, 3] <- x$scale
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
  spline_def = NULL
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
  noncentered_xi <- as.integer(spline_def$noncentered_xi %||% 0)

  c(
    list(
      use_alpha = as.integer(setup$has_icpt),
      use_beta = as.integer(setup$has_x),
      use_gamma = as.integer(setup$has_w),
      use_nu = as.integer(setup$has_w && spline_def$type == "rw2"),
      use_ar1 = as.integer(setup$has_ar1),
      use_zeta = as.integer(setup$has_zeta),
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
      likelihood = as.integer(!setup$prior_only),
      noncentered_xi = noncentered_xi,
      missing_idx = setup$missing_idx,
      n_missing = nrow(setup$missing_idx)
    ),
    unlist(
      lapply(
        c("sigma", "omega", "alpha", "beta", "kappa", "nu", "rho", "zeta"),
        \(i) prior_to_stan(priors[[i]], i)
      ),
      recursive = FALSE
    )
  )
}
# x: standata, d: descriptives
create_inits <- function(x, d, spline_def, marginalize_sigma = FALSE) {
  eta <- matrix(stats::rnorm(x$N * (x$J - 1), 0, 0.1), x$N, x$J - 1)
  sigma <- array(stats::runif(x$N, 0.5, 1.5) * d$sd_e)
  a <- beta <- rho <- kappa <- nu <- zeta <- array(0, 0)
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
      if (x$use_nu) {
        nu <- array(stats::rnorm(x$L, 0, 0.05 * sd_yx / d$sd_t))
      }
    }
  }
  if (x$use_ar1) {
    rho <- array(stats::runif(x$N, 0, 0.5))
  }
  if (x$use_zeta) {
    zeta <- array(stats::runif(x$N, 0.5, 1))
  }
  inits <- dplyr::lst(eta, sigma, a, beta, xi, kappa, nu, rho, zeta)
  # sigma is not a parameter of the marginalized model
  if (marginalize_sigma) {
    inits$sigma <- NULL
  }
  inits
}
