#' Bayesian R-squared value
#'
#' `bayes_R2` computes the Bayesian \eqn{R^2} measure of model fit for a 
#' Bayesian synthetic control model, while `loo_R2` computes a leave-one-out 
#' adjusted version of the same quantity.
#'
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return `data.frame` of posterior summary of (LOO-adjusted) R-squared values.
#' @references
#' Gelman A, Goodrich B, Gabry J, and Vehtari A (2019). R-squared for
#' Bayesian regression models. *The American Statistician*. 73(3), 307--309.
#' @rdname R2
#' @aliases loo_R2
#' @export loo_R2
#' @export
loo_R2.bscmfit <- function(object, probs = c(0.025, 0.975), ...) {
  stopifnot_(
    !is.null(object$data),
    "LOO R2 requires the original data. Refit the model with
    {.code save_data = TRUE}."
  )
  test_probs(probs)
  old_seed <- .Random.seed
  on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
  set.seed(get_stanfit(object)@stan_args[[1]]$seed)
  
  T_pre <- get_T_pre(object)
  unit <- get_unit(object)
  time <- get_time(object)
  times <- get_times(object)
  outcome <- get_outcome(object)
  treated <- get_treated(object)
  N <- get_N(object)
  mu_y <- posterior_epred(object)
  log_ratios <- -log_lik(object)
  cs <- cumsum(T_pre)
  idx1 <- stats::setNames(1L + c(0L, cs[-N]), treated)
  idx2 <- stats::setNames(cs, treated)
  r2 <- lapply(
    treated, \(id) {
      y <- object$data |>
        filter(.data[[unit]] == .env$id) |>
        pull(.data[[outcome]])
      y <- y[seq_len(T_pre[id])]
      idx <- idx1[id]:idx2[id]
      lr <- log_ratios[, idx, drop = FALSE]
      psis_object <- loo::psis(lr)
      mu_y_loo <- loo::E_loo(
        mu_y[, idx, drop = FALSE], psis_object, log_ratios = lr
      )$value
      e_loo <- mu_y_loo - y
      # Dirichlet weights for Bayesian bootstrap
      S <- nrow(mu_y)
      T_pre <- length(idx)
      w <- matrix(stats::rexp(S * T_pre), nrow = S, ncol = T_pre)
      w <- t(w / rowSums(w))
      ss_y <- (colSums(w * y^2) - colSums(w * y)^2)
      ss_e_loo <- (colSums(w * e_loo^2) - colSums(w * e_loo)^2)
      r2 <- 1 - ss_e_loo / ss_y
      r2[r2 < -1] <- -1
      r2[r2 > 1] <- 1
      posterior::draws_array(R2 = r2) |> 
        summarise_with_probs(probs)
    }
  ) |> stats::setNames(treated)
  bind_rows(r2, .id = unit)
}
#' @rdname R2
#' @aliases bayes_R2
#' @export bayes_R2
#' @export
bayes_R2.bscmfit <- function(object, probs = c(0.025, 0.975), ...) {
  test_probs(probs)
  
  d <- as_draws(object, "R2")
  N <- get_N(object)
  treated <- get_treated(object)
  unit <- get_unit(object)
  as_draws(object, "R2") |> 
    summarise_with_probs(probs) |>
    mutate("{unit}" := treated, .before = 1L) |> 
    mutate(variable = "R2")
}
