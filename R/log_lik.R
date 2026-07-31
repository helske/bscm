#' Posterior draws of pointwise log-likelihood
#'
#' Returns posterior draws of pointwise log-likelihoods of the treated units per
#' pre-treatment time point
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A matrix where rows correspond to posterior draws and columns to
#' observed (non-missing) time periods. Output columns are ordered by unit
#' then time: all observed pre-treatment time points for the first treated unit,
#' then the second, and so on.
#' @aliases log_lik
#' @export log_lik
#' @export
#' @examples
#' log_lik(fit_single_treated) |> head()
log_lik.bscmfit <- function(object, ...) {
  stopifnot_(
    !is.null(object$data),
    "Computing log-likelihood requires the original data. Refit the model
      with {.code save_data = TRUE}."
  )

  T_pre <- get_T_pre(object)
  T_total <- get_T_total(object)
  N <- get_N(object)
  n_obs <- sum(T_pre) - nrow(object$setup$missing_idx)
  mu_draws <- posterior_epred(object)
  sigma_draws <- as.matrix(get_stanfit(object), pars = "sigma")
  y <- get_stan_y(object)
  ll <- matrix(NA, nrow(mu_draws), n_obs)
  cn <- character(n_obs)
  idx <- 0L
  for (i in seq_len(N)) {
    mu_offset <- (i - 1L) * T_total
    sigma <- sigma_draws[, i]
    for (t in seq_len(T_pre[i])) {
      if (!is.na(y[t, i])) {
        idx <- idx + 1L
        ll[, idx] <- stats::dnorm(
          y[t, i],
          mu_draws[, mu_offset + t],
          sigma,
          log = TRUE
        )
        cn[idx] <- paste0("log_lik[", t, ",", i, "]")
      }
    }
  }
  colnames(ll) <- cn
  ll
}

loglik_at_times <- function(fit, T_start, T_end = T_start) {
  T_total <- get_T_total(fit)
  N <- get_N(fit)
  mu_draws <- posterior_epred(fit)
  sigma_draws <- as.matrix(get_stanfit(fit), pars = "sigma")
  y <- t(fit$standata$y)
  T_seq <- seq.int(T_start, T_end)
  ll <- matrix(0, nrow(mu_draws), length(T_seq))
  for (i in seq_len(N)) {
    mu_offset <- (i - 1L) * T_total
    sigma <- sigma_draws[, i]
    for (j in seq_along(T_seq)) {
      t <- T_seq[j]
      y_it <- y[t, i]
      if (!is.na(y_it)) {
        ll[, j] <- ll[, j] +
          stats::dnorm(
            y_it,
            mu_draws[, mu_offset + t],
            sigma,
            log = TRUE
          )
      }
    }
  }
  ll
}
