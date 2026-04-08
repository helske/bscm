#' Posterior Predictive Distribution of Bayesian Synthetic Control Model
#' 
#' Returns draws from the posterior predictive distribution of the Bayesian 
#' synthetic control. Note that this function does not simulate new 
#' realizations from this distribution, but rather returns the posterior draws 
#' computed during the model estimation.
#' 
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A matrix of posterior predictive draws of the synthetic control, 
#' where rows correspond to posterior draws and columns to time periods.
#' @aliases posterior_predict
#' @export posterior_predict
#' @export
posterior_predict.bscmfit <- function(object, ...) {
  as.matrix(get_stanfit(object), pars = "y_rep")
}
#' Posterior Draws of the Expected Predictive Distribution
#'
#' Returns draws from the posterior distribution of the expected value of
#' the Bayesian synthetic control model.
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A matrix where rows correspond to posterior draws and columns to
#'   time periods.
#' @aliases posterior_epred
#' @export posterior_epred
#' @export
posterior_epred.bscmfit <- function(object, ...) {
  as.matrix(get_stanfit(object), pars = "y_mean")
}
#' Posterior Draws of the Linear Predictor
#'
#' Returns draws from the posterior distribution of the linear predictor of
#' the Bayesian synthetic control model. Since the model uses a Gaussian
#' likelihood with identity link, this is equivalent to
#' [posterior_epred.bscmfit()].
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param transform \[`logical(1)`]\cr Ignored.
#' @param ... Ignored.
#' @return A matrix where rows correspond to posterior draws and columns to
#' time periods.
#' @aliases posterior_linpred
#' @export posterior_linpred
#' @export
posterior_linpred.bscmfit <- function(object, transform = FALSE, ...) {
  as.matrix(get_stanfit(object), pars = "y_mean")
}
#' Posterior draws of pointwise log-likelihood
#'
#' Returns draws of log-likelihood values of the treated unit in the 
#' time points corresponding to the pre-treatment period.
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A matrix where rows correspond to posterior draws and columns to
#' time periods.
#' @aliases log_lik
#' @export log_lik
#' @export
log_lik.bscmfit <- function(object, ...) {
  stopifnot_(
    !is.null(object$data),
    "Computing log-likelihood requires the original data. Refit the model
    with {.code save_data = TRUE}."
  )
  T_total <- get_T_total(object)
  T_pre <- get_T_pre(object)
  N <- get_N(object)
  treated <- get_treated(object)
  unit <- get_unit(object)
  outcome <- get_outcome(object)
  mu_draws <- posterior_epred(object)
  sigma_draws <- as.matrix(get_stanfit(object), pars = "sigma")
  n_ll <- sum(T_pre)
  ll <- matrix(NA_real_, nrow(mu_draws), n_ll)
  idx <- 0L
  for (i in seq_len(N)) {
    id <- treated[i]
    y_obs <- object$data |>
      filter(.data[[unit]] == .env$id) |>
      pull(.data[[outcome]])
    y_pre <- y_obs[seq_len(T_pre[id])]
    mu_offset <- (i - 1L) * T_total
    s <- sigma_draws[, i]
    for (t in seq_along(y_pre)) {
      idx <- idx + 1L
      ll[, idx] <- stats::dnorm(
        y_pre[t], mu_draws[, mu_offset + t], s, log = TRUE
      )
    }
  }
  colnames(ll) <- paste0("log_lik[", seq_len(n_ll), "]")
  ll
}
