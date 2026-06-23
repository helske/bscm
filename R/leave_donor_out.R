#' @export
#' @rdname leave_donor_out
leave_donor_out <- function(x, ...) {
  UseMethod("leave_donor_out", x)
}
#' Leave-out donor sensitivity of a Bayesian synthetic control model
#'
#' `leave_donor_out()` re-estimates the original model after omitting donors
#' from the donor pool. By default, it removes one donor at a time. Setting
#' `cumulative = TRUE` instead removes donors cumulatively following the chosen
#' order, which can be used to assess how sensitive the results are to the most
#' influential donors.
#'
#' The donor order can always be supplied explicitly. For one-donor leave-out
#' runs (`cumulative = FALSE`), that order controls the order of the returned
#' runs. For cumulative leave-out runs (`cumulative = TRUE`), donors are removed
#' from the ordered pool one by one. If no order is supplied, one-donor
#' leave-out uses the original donor order, while cumulative leave-out orders
#' donors by posterior mean donor weight in descending order.
#'
#' For models with multiple treated units, automatic donor ranking is based on
#' the average posterior mean donor weight across treated units.
#'
#' @param x \[`bscmfit`]\cr The output returned by the [bscm()].
#' @param order \[`character()` or `NULL`]\cr Donor order used for leave-out
#'   runs. Either a vector of donor names, or `"descending"` / `"ascending"` to
#'   rank donors by posterior means of donor weights. Default is `"descending"`.
#'   Use `NULL` to use the original donor order.
#' @param cumulative \[`logical(1)`]\cr If `FALSE` (the default), omit one donor
#'   at a time. If `TRUE`, omit donors cumulatively following `order`.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries of the
#' treatment effects and RMSE estimates. Default is `c(0.025, 0.975)`.
#' @param ... Additional arguments passed on to [bscm()].
#' @return An object of class `bscm_ldo` with data frames `effect`,
#'   `rmse`, `weights`, and `diagnostics`, and a `metadata` list. The data
#'   frames contain posterior summaries for the original fit (`step = 0`) and
#'   for each leave-out run, identified by `step`, `n_removed`, and
#'   `last_removed`. The `metadata` list contains the omitted donors, whether
#'   donor omission was cumulative, summary probabilities, and model metadata
#'   needed for plotting. The result can be visualized with [plot_weights()] and
#'   [plot_effects()].
#' @rdname leave_donor_out
#' @aliases leave_donor_out
#' @export
leave_donor_out.bscmfit <- function(
  x,
  order = "descending",
  cumulative = FALSE,
  probs = c(0.025, 0.5, 0.975),
  ...
) {
  probs <- sort_probs(probs)
  stopifnot_(
    checkmate::test_flag(cumulative),
    "Argument {.arg cumulative} must be either TRUE or FALSE."
  )
  stopifnot_(
    !is.null(x$data),
    "The model fit {.arg x} does not contain the original data. You probably
    used {.fun bscm} with {.arg save_data = FALSE}?"
  )
  stopifnot_(
    get_J(x) >= 2L,
    "Donor exclusion requires at least two donors in the donor pool."
  )
  donor_order <- order_donors(x, order)
  if (!cumulative) {
    removed_donors <- stats::setNames(as.list(donor_order), donor_order)
  } else {
    removed_donors <- donor_order[-length(donor_order)]
    removed_donors <- stats::setNames(
      lapply(seq_along(removed_donors), \(i) utils::head(removed_donors, i)),
      removed_donors
    )
  }
  unit <- get_unit(x)
  effects <- weights <- rmses <- diagnostics <-
    stats::setNames(
      vector("list", length(removed_donors) + 1L),
      c("none", names(removed_donors))
    )
  effects[[1]] <- treatment_effect(x, probs = probs)
  rmses[[1]] <- rmse(x, probs = probs)
  weights[[1]] <- donor_weights(x, probs = probs)
  diagnostics[[1]] <- check_mcmc_diagnostics(x, warn = FALSE)

  priors <- get_priors(x)
  p <- progressr::progressor(along = removed_donors)
  for (i in seq_along(removed_donors)) {
    p(paste0("Removing donor ", donor_order[i + 1]))
    d <- x$data |>
      dplyr::filter(!(.data[[unit]] %in% .env$removed_donors[[i]]))
    fit <- stats::update(
      x,
      data = d,
      refresh = 0,
      mcmc_diagnostics = FALSE,
      priors = priors,
      ...
    )
    effects[[i + 1]] <- treatment_effect(fit, probs = probs)
    rmses[[i + 1]] <- rmse(fit, probs = probs)
    weights[[i + 1]] <- donor_weights(fit, probs = probs)
    diagnostics[[i + 1]] <- check_mcmc_diagnostics(fit, warn = FALSE)
  }
  issues <- vapply(
    diagnostics[-1L],
    \(x) x$has_issues,
    logical(1)
  )
  warnifnot_(
    all(!issues),
    "Some of the runs resulted in MCMC diagnostic warnings. Check
    the `diagnostics` element of the output list for details."
  )
  out <- list(
    effect = dplyr::bind_rows(effects, .id = "removed_donor"),
    rmse = dplyr::bind_rows(rmses, .id = "removed_donor"),
    weights = dplyr::bind_rows(weights, .id = "removed_donor"),
    diagnostics = dplyr::bind_rows(
      lapply(diagnostics, diags2df),
      .id = "removed_donor"
    ),
    metadata = list(
      removed_donors = removed_donors,
      cumulative = cumulative,
      donor_order = donor_order,
      order = order,
      probs = probs,
      time = x$setup$time,
      unit = unit,
      treated = x$setup$treated
    )
  )
  class(out) <- "bscm_ldo"
  out
}

diags2df <- function(x) {
  data.frame(
    has_issues = x$has_issues,
    n_divergences = x$n_divergences,
    n_max_treedepth = x$n_max_treedepth,
    n_low_bfmi = x$n_low_bfmi,
    messages = I(list(x$messages)),
    rhat_and_ess = I(list(x$rhat_and_ess))
  )
}
