#' Posterior draws of the Bayesian R-squared value
#'
#' `bayes_R2` computes the Bayesian \eqn{R^2} measure of model fit for a 
#' Bayesian synthetic control model, while `loo_R2` computes a leave-one-out 
#' adjusted version of the same quantity.
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param summary \['logical(1)]\cr If `TRUE` (the default) returns posterior 
#'   summary statistics, if `FALSE`, returns posterior draws.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return Either a vector of posterior draws or `data.frame` of posterior 
#'   summary of (LOO-adjusted) R-squared values.
#' @references
#' Gelman A, Goodrich B, Gabry J, and Vehtari A (2019). R-squared for
#' Bayesian regression models. *The American Statistician*. 73(3), 307--309.
#' @aliases loo_R2
#' @rdname bayes_R2
#' @export
#' @export loo_R2
loo_R2.bscmfit <- function(object, summary = TRUE, 
                           probs = c(0.025,0.975), ...) {
  stopifnot_(
    !is.null(object$data),
    "LOO R-squared requires the original data. Refit the model with
    {.code save_data = TRUE}."
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
  T_pre <- get_T_pre(object)
  y <- object$data |>
    filter(.data[[get_unit(object)]] == get_treated(object)) |>
    pull(.data[[get_outcome(object)]])
  y <- y[seq_len(T_pre)]
  mu_y <- posterior_epred(object)[, seq_len(T_pre), drop = FALSE]
  log_ratios <- -log_lik(object)
  psis_object <- loo::psis(log_ratios)
  mu_y_loo <- loo::E_loo(mu_y, psis_object, log_ratios = log_ratios)$value
  e_loo <- mu_y_loo - y
  
  old_seed <- .Random.seed
  on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
  set.seed(object$stanfit@stan_args[[1]]$seed)
  
  # Dirichlet weights for Bayesian bootstrap
  S <- nrow(mu_y)
  w <- matrix(stats::rexp(S * T_pre), nrow = S, ncol = T_pre)
  w <- t(w / rowSums(w))
  
  ss_y <- (colSums(w * y^2) - colSums(w * y)^2)
  ss_e_loo <- (colSums(w * e_loo^2) - colSums(w * e_loo)^2)
  
  r2 <- 1 - ss_e_loo / ss_y
  r2[r2 < -1] <- -1
  r2[r2 > 1] <- 1
  if (summary) {
    posterior::draws_array(R2 = r2) |> 
      summarize_with_probs(probs)
  } else {
    r2
  }
}
#' @aliases bayes_R2
#' @rdname bayes_R2
#' @export
#' @export bayes_R2
bayes_R2.bscmfit <- function(object, summary = TRUE, 
                             probs = c(0.025,0.975), ...) {
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
  r2 <- c(as.matrix(object$stanfit, "R2"))
  if (summary) {
    posterior::draws_array(R2 = r2) |> 
      summarize_with_probs(probs)
  } else {
    r2
  }
}