#' @export
#' @rdname covariate_distance
covariate_distance <- function(x, ...) {
  UseMethod("covariate_distance", x)
}
#' Synthetic covariate distance
#'
#' For models with covariates, returns and optionally visualizes the synthetic 
#' covariate distances \eqn{\sqrt{\delta_{k,t}}}, 
#' where
#' \deqn{\delta_{k,t} = (x_{k,0,t} - \sum_{j=1}^J \omega_j x_{k, j, t})^2 / Var(x_k),}
#' \eqn{t=1,\ldots,T}, \eqn{x_{k,0,t}} is the value of $k$th covariate of 
#' treated unit at time t, and similarly for donors \eqn{j=1,\ldots,J}.
#' 
#' @export
#' @rdname covariate_distance
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param plot \[`logical(1)`]\cr If `TRUE` (the default), plots the posterior 
#' mean and interval of the synthetic covariate distances over time.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#' Default is `c(0.025, 0.975)`. If length of `probs` less than 2, no posterior 
#' intervals are drawn, and if length of `probs` is larger than two, the most 
#' extreme values are used for the posterior intervals.
#' @param ... Optional arguments passed to [ggplot2::facet_wrap()].
#' @return A `data.frame` of posterior summaries of synthetic covariate
#' distances.
covariate_distance.bscmfit <- function(x, plot = TRUE, probs = c(0.025, 0.975),
                                       ...) {
  variable <- NULL # to avoid NSE warnings
  stopifnot_(
    checkmate::test_flag(plot),
    "Argument {.arg plot} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_numeric(
      probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values between
     0 and 1."
  )
  
  stopifnot_(
    has_predictors(x),
    "The model does not contain any predictors."
  )
  X <- get_Xs(x)
  omega <- as_draws_rvars(as_draws(x, "omega"))$omega
  delta <- delta_mean <- stats::setNames(
    vector("list", ncol(X$X_y)),
    colnames(X$X_y) 
  )
  time <- get_time(x)
  times <- get_times(x)
  for (k in seq_along(delta)) {
    dx <- (X$X_y[, k] - c(omega %*% X$X_z[, , k]))^2
    delta[[k]] <- dx |> 
      summarise_draws(
        mean, sd, ~ quantile2(.x, probs = probs), 
        default_convergence_measures()
      ) |> 
      mutate("{time}" := .env$times, .before = variable) |> 
      select(-variable)
    delta_mean[[k]] <- rvar_mean(dx) |>  
      summarise_draws(mean)
    
  }
  delta <- bind_rows(delta, .id = "Covariate")
  delta_mean <- bind_rows(delta_mean, .id = "Covariate")
  if (plot) {
    lwr <- paste0("q", 100 * min(probs))
    upr <- paste0("q", 100 * max(probs))
    if (lwr != upr) {
      ribbon <- geom_ribbon(
        aes(ymin = .data[[lwr]], ymax = .data[[upr]]), alpha = 0.25
      )
    } else {
      ribbon <- NULL
    }
    p <- delta |> 
      ggplot(aes(.data[[time]], mean)) +
      ribbon +
      geom_hline(
        data = delta_mean, aes(yintercept = mean), 
        linetype = "dashed", colour = "grey50"
      ) +
      geom_line() +
      labs(x = time, y = "Synthetic covariate distance") +
      facet_wrap( ~ Covariate, ...) +
      theme_bw()
    print(p)
  }
  delta
}

get_Xs <- function(x) {
  X <- stats::model.matrix(stats::formula(x), data = x$data)
  if (has_intercept(x)) X <- X[, -1L, drop = FALSE]
  K <- ncol(X)
  T_total <- get_T_total(x)
  T_pre <- get_T_pre(x)
  J <- length(get_donors(x))
  X <- simplify2array(
    lapply(seq_len(K), \(k) matrix(X[, k], J + 1, T_total, TRUE))
  )
  sd_x <- apply(
    apply(X[, 1:T_pre, , drop = FALSE], c(1, 3), sd), 
    2, stats::median
  )
  for (k in seq_len(K)) {
    X[, , k] <- X / sd_x[k]
  }
  treated_idx <- which(unique(x$data[[get_unit(x)]]) == get_treated(x))
  X_z <- X[-treated_idx, , , drop = FALSE]
  X_y <- X[treated_idx, , , drop = FALSE]
  dim(X_y) <- c(T_total, K)
  list(X_y = X_y, X_z = X_z)
}
