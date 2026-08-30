#' Bayesian R-squared value
#'
#' `bayes_R2` computes the Bayesian \eqn{R^2} measure of model fit for a
#' Bayesian synthetic control model, while `loo_R2` computes a leave-one-out
#' adjusted version of the same quantity.
#'
#' @inheritParams bscm_postprocessing
#' @param fixed_seed \[logical(1)]\cr If `TRUE` (the default), fixes the seed
#'   of random number generator (RNG) so that `loo_R2`, which uses Bayesian
#'   bootstrap, returns identical results in repeated calls for the same model
#'   object. On exit, the state of the RNG is restored to the original state.
#' @return A `tibble` of posterior summary of (LOO-adjusted) R-squared values.
#' @references
#' Gelman A, Goodrich B, Gabry J, and Vehtari A (2019). R-squared for
#' Bayesian regression models. *The American Statistician*. 73(3), 307--309,
#' <doi:10.1080/00031305.2018.1549100>.
#' @rdname R2
#' @aliases loo_R2
#' @export loo_R2
#' @export
#' @examples
#' bayes_R2(fit_single_treated)
loo_R2.bscmfit <- function(
  object,
  summary = TRUE,
  probs = c(0.025, 0.975),
  fixed_seed = TRUE,
  ...
) {
  check_flag(summary, "summary")
  check_flag(fixed_seed, "fixed_seed")
  probs <- sort_probs(probs)
  check_has_data(object, "object")
  if (fixed_seed) {
    old_seed <- .Random.seed
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
    set.seed(get_stanfit(object)@stan_args[[1]]$seed)
  }
  T_pre <- get_T_pre(object)
  unit <- get_unit(object)
  outcome <- get_outcome(object)
  treated <- get_treated(object)
  N <- get_N(object)
  mu_y <- posterior_epred(object)
  log_ratios <- -log_lik(object)
  cs <- cumsum(T_pre)
  idx1 <- stats::setNames(1L + c(0L, cs[-N]), treated)
  idx2 <- stats::setNames(cs, treated)
  d <- lapply(
    treated,
    \(id) {
      y <- object$data |>
        dplyr::filter(.data[[unit]] == .env$id) |>
        dplyr::pull(.data[[outcome]])
      y <- y[seq_len(T_pre[id])]
      idx <- idx1[id]:idx2[id]
      lr <- log_ratios[, idx, drop = FALSE]
      psis_object <- loo::psis(lr)
      mu_y_loo <- loo::E_loo(
        mu_y[, idx, drop = FALSE],
        psis_object,
        log_ratios = lr
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
      dplyr::tibble(
        "{unit}" := id,
        R2 = posterior::draws_rvars(r2 = r2)$r2
      )
    }
  ) |>
    dplyr::bind_rows()
  if (summary) {
    d <- summarise_column(d, "R2", probs)
  }
  d
}
#' @rdname R2
#' @aliases bayes_R2
#' @export bayes_R2
#' @export
bayes_R2.bscmfit <- function(
  object,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ...
) {
  check_flag(summary, "summary")
  probs <- sort_probs(probs)
  N <- get_N(object)
  treated <- get_treated(object)
  unit <- get_unit(object)
  T_pre <- get_T_pre(object)

  y_mean <- gq_to_rvar(object, object$y_mean)
  sigma <- rvars_of(object, "sigma")

  r2 <- vector("list", N)
  for (i in seq_len(N)) {
    var_fit <- posterior::rvar_var(y_mean[seq_len(T_pre[treated[i]]), i])
    r2[[i]] <- var_fit / (var_fit + sigma[i]^2)
  }
  d <- dplyr::tibble(
    "{unit}" := treated,
    R2 = do.call(c, r2)
  )
  if (summary) {
    d <- summarise_column(d, "R2", probs)
  }
  d
}
