#' @export
#' @rdname covariate_adjustment
covariate_adjustment <- function(x, ...) {
  UseMethod("covariate_adjustment", x)
}
#' Covariate adjustments for Bayesian synthetic control model
#'
#' Computes the posterior distribution of the covariate adjustment to the
#' synthetic control for each predictor, defined as the
#' regression coefficient multiplied by the difference between the predictor
#' value of the treated unit and the corresponding weighted predictor value
#' of the donor units.
#'
#' For predictor \eqn{k}, treated unit \eqn{i}, and time point \eqn{t}, the
#' covariate adjustment is
#' \deqn{
#'   b_{tk}
#'   \left(
#'     X_{y,itk} - \sum_j \omega_{ji} X_{z,jtk}
#'   \right),
#' }
#' where \eqn{b_{tk} = \beta_k} for predictors with a fixed coefficient and
#' \eqn{b_{tk} = \beta_k + \gamma_{tk}} for predictors with a time-varying
#' coefficient.
#'
#' If there are multiple treated units, the covariate adjustments are averaged
#' over treated units within each posterior draw before computing posterior
#' summaries.
#'
#' @inheritParams plot_coefs
#' @param x \[`bscmfit`]\cr An object of class `bscmfit`.
#' @param plot \[`logical(1)`]\cr If `TRUE` (the default), plot the posterior
#'   covariate adjustments.
#' @param ... Ignored.
#' @return A data frame containing the posterior mean and requested posterior
#'   quantiles of the covariate adjustment for each predictor and time point.
#'   If `plot = TRUE`, the corresponding plot is printed as a side effect.
#' @rdname covariate_adjustment
#' @export
covariate_adjustment.bscmfit <- function(
    x,
    plot = TRUE,
    probs = c(0.025, 0.975),
    alpha = 0.5,
    scales = "free_y",
    ...
) {
  stopifnot_(
    has_predictors(x),
    "The model does not contain any predictors."
  )
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
      len = 2L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector of length 2 with
    values between 0 and 1."
  )
  stopifnot_(
    checkmate::test_number(
      alpha,
      lower = 0.0,
      upper = 1.0
    ),
    "Argument {.arg alpha} must be a single {.cls numeric} value between 0
    and 1."
  )
  
  probs <- sort(probs)
  X <- get_Xs(x)
  N <- get_N(x)
  time_var <- get_time(x)
  times <- get_times(x)
  predictors <- get_predictors(x)
  
  adjustment <- Reduce(
    `+`,
    lapply(seq_len(N), \(i) covariate_adjustment_unit(x, i, X))
  ) / N
  
  d <- lapply(
    seq_along(predictors),
    \(k) {
      adjustment[, k] |>
        summarise_with_probs(probs, for_plots = TRUE) |>
        dplyr::mutate(
          "{time_var}" := .env$times,
          variable = .env$predictors[k],
          .before = 1L
        )
    }
  ) |>
    dplyr::bind_rows()
  
  if (plot) {
    qs <- paste0("q", 100 * probs)
    
    p <- d |>
      ggplot(aes(.data[[time_var]], mean)) +
      geom_ribbon(
        aes(ymin = .data[[qs[1]]], ymax = .data[[qs[2]]]),
        fill = "#77AADD",
        alpha = alpha
      ) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_line(colour = "#0C7BDC") +
      scale_x_continuous(limits = range(d[[time_var]])) +
      facet_wrap(~variable, scales = scales, ...) +
      labs(x = time_var, y = "Covariate adjustment") +
      theme_bw()
    
    print(p)
  }
  d
}
#' @noRd
covariate_adjustment_unit <- function(x, unit, X) {
  T_total <- get_T_total(x)
  J <- get_J(x)
  K <- length(get_predictors(x))
  predictors <- get_predictors(x)
  
  pars <- paste0("omega[", seq_len(J), ",", unit, "]")
  omega <- posterior::rvar(
    as_draws(x, pars),
    with_chains = TRUE
  )
  
  d_beta <- coef(x, type = "beta", summary = FALSE)$beta |>
    dplyr::mutate(
      parameter = sub("^beta_", "", .data$parameter)
    )
  
  if (has_tv_coefs(x)) {
    d_gamma <- coef(x, type = "gamma", summary = FALSE)$gamma |>
      dplyr::mutate(
        parameter = sub("^gamma_", "", .data$parameter)
      )
    
    d_beta <- d_beta |>
      dplyr::left_join(d_gamma, by = "parameter") |>
      dplyr::mutate(
        beta = dplyr::coalesce(
          .data$beta + .data$gamma,
          .data$beta
        )
      ) |>
      dplyr::select(-"gamma")
  }
  
  X_y <- X$X_y[unit, , , drop = FALSE]
  dim(X_y) <- c(T_total, K)
  
  adjustment <- lapply(
    seq_len(K),
    \(k) {
      delta <- X_y[, k] - c(omega %*% X$X_z[, , k])
      b <- d_beta$beta[d_beta$parameter == predictors[k]]
      b * delta
    }
  )
  
  do.call(cbind, adjustment)
}
