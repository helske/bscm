#' Approximate leave-one-out (LOO) cross-validation for Bayesian synthetic
#' control models
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param r_eff \[`logical(1)`]\cr If `TRUE` (the default), [loo::loo()]
#' computes more accurate Monte Carlo error estimates at the cost of
#' increased computation time.
#' @param reloo \[`logical(1)`]\cr If `TRUE`, observations whose Pareto
#' `k` diagnostic exceeds `k_threshold` are handled by exact
#' leave-one-out refits rather than PSIS approximation. Requires
#' `save_data = TRUE` in the original fit. Default is `FALSE`.
#' @param k_threshold \[`numeric(1)`]\cr Threshold for the Pareto `k`
#' diagnostic that triggers an exact refit when `reloo = TRUE`.
#' Default is `0.7`.
#' @param ... Additional arguments to [loo::loo()].
#' @return An output from [loo::loo()].
#' @references
#' Vehtari A, Gelman A, and Gabry J (2017).
#' Practical Bayesian model evaluation using leave-one-out cross-validation and
#' WAIC. *Statistics and Computing*. 27(5), 1413--1432,
#' <doi:10.1007/s11222-016-9696-4>.
#' @aliases loo
#' @export loo
#' @export
#' @examples
#' loo(fit_single_treated)
loo.bscmfit <- function(x, 
                        r_eff = TRUE, 
                        reloo = FALSE, 
                        k_threshold = 0.7, 
                        ...) 
{
  
  ll <- log_lik(x)
  nc <- nchains(x)
  cid <- rep(seq_len(nc), each = ndraws(x) / nc)
  r_eff_val <- if (r_eff) loo::relative_eff(exp(ll), chain_id = cid) else 1
  loo_out <- loo::loo(ll, r_eff = r_eff_val, ...)
  if (!reloo) {
    return(loo_out)
  }
  bad_k <- loo::pareto_k_ids(loo_out, threshold = k_threshold)
  if (length(bad_k) == 0L) {
    return(loo_out)
  }
  T_total <- get_T_total(x)
  standata <- get_standata(x)
  standata$sample_y_rep <- FALSE
  ti <- matrix(
    as.integer(
      unlist(
        strsplit(
          sub("^log_lik\\[|\\]$", "", colnames(ll)[bad_k]),
          ",",
          fixed = TRUE
        )
      )
    ),
    2,
    length(bad_k)
  )
  for (k in seq_along(bad_k)) {
    obs <- bad_k[k]
    t <- ti[1, k]
    i <- ti[2, k]
    standata_mod <- standata
    standata_mod$missing_idx <- rbind(
      standata$missing_idx,
      c(unit = i, time = t)
    )
    standata_mod$n_missing <- standata$n_missing + 1L
    refit <- refit_bscm(x, standata_mod)
    mu_offset <- (i - 1L) * T_total
    mu_draws <- posterior_epred(refit)
    sigma_draws <- as.matrix(get_stanfit(refit), pars = "sigma")
    ll_exact <- stats::dnorm(
      standata$Y[i, t],
      mu_draws[, mu_offset + t],
      sigma_draws[, i],
      log = TRUE
    )
    elpd_exact <- log_mean_exp(ll_exact)
    lpd_i <- log_mean_exp(ll[, obs])
    loo_out$pointwise[obs, "elpd_loo"] <- elpd_exact
    loo_out$pointwise[obs, "looic"] <- -2 * elpd_exact
    loo_out$pointwise[obs, "p_loo"] <- lpd_i - elpd_exact
    loo_out$diagnostics$pareto_k[obs] <- 0
  }
  pw <- loo_out$pointwise[, c("elpd_loo", "p_loo", "looic"), drop = FALSE]
  loo_out$estimates[, "Estimate"] <- colSums(pw)
  loo_out$estimates[, "SE"] <-
    sqrt(nrow(pw) * apply(pw, 2, stats::var))
  loo_out
}
