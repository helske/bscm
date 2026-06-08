#' Bayesian R-squared value
#'
#' `bayes_R2` computes the Bayesian \eqn{R^2} measure of model fit for a
#' Bayesian synthetic control model, while `loo_R2` computes a leave-one-out
#' adjusted version of the same quantity.
#'
#' @inheritParams rmse.bscmfit
#' @param object \[`bscmfit`]\cr The model fit object.
#' @param fixed_seed \[logical(1)]\cr If `TRUE` (the default), fixes the seed
#'   of random number generator (RNG) so that `loo_R2`, which uses Bayesian
#'   bootstrap, returns identical results in repeated calls for the same model
#'   object. On exit, the state of the RNG is restored to the original state.
#' @param ... Ignored.
#' @return `data.frame` of posterior summary of (LOO-adjusted) R-squared values.
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
  stopifnot_(
    !is.null(object$data),
    "LOO R2 requires the original data. Refit the model with
    {.code save_data = TRUE}."
  )
  test_summary(summary)
  probs <- sort_probs(probs)
  stopifnot_(
    checkmate::test_flag(fixed_seed),
    "Argument {.arg fixed_seed} should be a single logical value."
  )
  if (fixed_seed) {
    old_seed <- .Random.seed
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
    set.seed(get_stanfit(object)@stan_args[[1]]$seed)
  }
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
      tibble(
        "{unit}" := id,
        R2 = posterior::draws_rvars(r2 = r2)$r2
      )
    }
  ) |>
    dplyr::bind_rows()
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(.data$R2, probs)) |>
      dplyr::select(-"R2", -"variable")
  }
  if (length(treated) == 1) {
    d <- d |> dplyr::select(-dplyr::all_of(unit))
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
  test_summary(summary)
  probs <- sort_probs(probs)
  N <- get_N(object)
  treated <- get_treated(object)
  unit <- get_unit(object)
  T_pre <- get_T_pre(object)

  y_mean <- posterior::as_draws_rvars(as_draws(object, "y_mean"))$y_mean
  sigma <- posterior::as_draws_rvars(as_draws(object, "sigma"))$sigma

  r2 <- vector("list", N)
  for (i in seq_len(N)) {
    var_fit <- posterior::rvar_var(y_mean[seq_len(T_pre[treated[i]]), i])
    r2[[i]] <- var_fit / (var_fit + sigma[i]^2)
  }
  d <- tibble(
    "{unit}" := treated,
    R2 = do.call(c, r2)
  )
  if (summary) {
    d <- d |>
      dplyr::mutate(summarise_with_probs(.data$R2, probs)) |>
      dplyr::select(-"R2", -"variable")
  }
  if (N == 1) {
    d <- d |> dplyr::select(-dplyr::all_of(unit))
  }
  d
}
