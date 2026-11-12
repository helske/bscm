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
    if (inherits(x$concentration, "bscm_prior")) {
      x$length <- 0L
      x$concentration <- NULL
    }
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
  names(priors) <- stan_parameter_names[names(priors)]
  if (inherits(priors$omega$concentration, "bscm_prior")) {
    priors$zeta <- priors$omega$concentration
  }
  
  standata <- c(
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
  tilt <- priors$omega$tilt
  if (!is.null(tilt)) {
    standata <- c(
      standata,
      create_tilt_data(tilt, Z, min(setup$T_pre), setup$donors)
    )
  }
  standata
}

#' Create the Stan data for the tilt of the prior of the donor weights
#'
#' With the tilt, the donor weights of treated unit `i` are
#' `softmax(V phi_i + ...)` with `phi_i ~ normal(0, scale)`. When the tilt is
#' defined by the number of directions `R`, the columns of `V` are the first
#' `R` right singular vectors of the pre-treatment donor outcomes after
#' removing the donor means and the common time effects, standardized to unit
#' standard deviation. A user-supplied `V` is used as is, with rows matched to
#' the donors by the row names.
#'
#' @param tilt \[`list()`]\cr Tilt specification from `check_tilt()`.
#' @param Z \[`matrix()`]\cr Donor outcomes, with time points on rows.
#' @param T_pre \[`integer()`]\cr Numbers of pre-treatment time points of the
#'   treated units.
#' @param donors \[`character()`]\cr Donor ids, in the order of the columns
#'   of `Z`.
#' @param call The calling environment reported in the error message.
#' @return A list with elements `R`, `V`, and `sd_phi`.
#' @noRd
create_tilt_data <- function(
    tilt,
    Z,
    T_pre,
    donors,
    call = rlang::caller_env()
) {
  J <- ncol(Z)
  T0 <- min(T_pre)
  if (is.null(tilt$V)) {
    R <- tilt$R
    max_R <- min(T0, J) - 1L
    stopifnot_(
      R <= max_R,
      c(
        "Element {.val R} of {.arg tilt} must be at most {max_R}.",
        i = "The tilt directions are computed from {T0} pre-treatment time
        point{?s} of {J} donors, so there are at most {max_R} of them.",
        i = "The default prior of {.val weights} has a tilt with {.code R = 3}.
        Define the prior of {.val weights} with smaller {.val R} or without a 
        tilt."
      ),
      call = call
    )
    Z_pre <- Z[seq_len(T0), , drop = FALSE]
    Z_pre <- sweep(Z_pre, 2, colMeans(Z_pre))
    Z_pre <- Z_pre - rowMeans(Z_pre)
    V <- unclass(scale(svd(Z_pre, nu = 0, nv = R)$v))
  } else {
    missing_donors <- setdiff(donors, rownames(tilt$V))
    stopifnot_(
      length(missing_donors) == 0L,
      c(
        "Element {.val V} of {.arg tilt} has no rows for
        {cli::qty(missing_donors)}donor{?s} {.val {missing_donors}}.",
        i = "The rows of {.val V} are matched to the donors by the row names."
      ),
      call = call
    )
    V <- unname(tilt$V[donors, , drop = FALSE])
  }
  list(R = ncol(V), V = V, sd_phi = tilt$scale)
}

# x: standata, d: descriptives
create_inits <- function(x, d, spline_def, marginalize_sigma = FALSE) {
  a <- beta <- rho <- kappa <- nu <- zeta <- array(0, 0)
  eta <- matrix(stats::rnorm(x$N * (x$J - 1), 0, 0.1), x$N, x$J - 1)
  sigma <- array(
    pmax(abs(stats::rnorm(x$N, 0, 0.5 * d$sd_e)), 0.01)
  )
  xi <- matrix(0, 0, 0)
  if (x$use_alpha) {
    a <- array(stats::rnorm(x$N, d$mean_y, 0.25 * d$sd_e))
  }
  if (x$use_beta) {
    sd_yx <- d$sd_y / d$sd_x
    beta <- array(stats::rnorm(x$K, 0, 0.1 * sd_yx))
    if (x$use_gamma) {
      xi <- matrix(rnorm(x$D * x$L, 0, 0.01), x$D, x$L)
      kappa <- array(
        pmax(abs(stats::rnorm(x$L, 0, 0.05 * sd_yx[x$tv_idx])), 0.01)
      )
      if (x$use_nu) {
        c_nu <- 2 * sqrt(2) / (max(x$T_pre) - 1)
        nu <- array(stats::rnorm(x$L, 0, 0.1 * c_nu * sd_yx[x$tv_idx]))
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
  # tilt coefficients exist only in the model with the tilt
  if (!is.null(x$V)) {
    inits$phi <- matrix(stats::rnorm(x$N * x$R, 0, 0.1), x$N, x$R)
  }
  inits
}
