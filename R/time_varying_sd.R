#' @export
#' @rdname time_varying_sd
time_varying_sd <- function(x, ...) {
  UseMethod("time_varying_sd", x)
}
#' Standard deviation induced by the time-varying covariate imbalance effects
#'
#' For models estimate with argument `time_varying_effects = TRUE`, 
#' returns and optionally visualizes the standard 
#' deviation due to the time-varying covariate imbalance effects, i.e., 
#' \eqn{\tau\sqrt{\delta_t}}, 
#' \deqn{\delta_t = \sum_{k=1}^K(x_{k,0,t} - \sum_{j=1}^J \omega_j x_{k, j, t})^2 / Var(x_k),}
#' \eqn{t=1,\ldots,T}, \eqn{x_{k,0,t}} is the value of _k_th covariate of 
#' treated unit at time t, and similarly for donors \eqn{j=1,\ldots,J}.
#' 
#' @export
#' @rdname time_varying_sd
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param plot \[`logical(1)`]\cr If `TRUE` (the default), plots the posterior 
#' mean and interval of the (scaled) synthetic covariate distances over time.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#' Default is `c(0.025, 0.975)`. If length of `probs` less than 2, no posterior 
#' intervals are drawn, and if length of `probs` is larger than two, the most 
#' extreme values are used for the posterior intervals.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the standard deviation of 
#' the time-varying covariate imbalance effects.
time_varying_sd.bscmfit <- function(x, plot = TRUE, probs = c(0.025, 0.975), 
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
    x$setup$time_varying_effects,
    "The model was not estimated with {.arg time_varying_effects = TRUE}."
  )
  delta <- sqrt(as_draws_rvars(as_draws(x, "delta"))$delta)
  delta <- delta * as_draws_rvars(as_draws(x, "tau"))$tau
  time <- get_time(x)
  times <- get_times(x)
  delta <- delta |> summarise_draws(
    mean, sd, 
    ~ quantile2(.x, probs = probs), 
    default_convergence_measures()
  ) |> 
    mutate("{time}" := .env$times, .before = variable) |> 
    select(-variable)
  
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
    sigma <- mean(as_draws(x, "sigma"))
    p <- delta |> 
      ggplot(aes(.data[[time]], mean)) +
      ribbon +
      geom_hline(
        aes(yintercept = sigma, linetype = "sigma")
      ) +
      geom_line(aes(linetype = "taudelta")) +
      scale_linetype_manual(
        name = NULL,
        values = c(
          taudelta = "solid",
          sigma = "dashed"
        ),
        labels = c(
          taudelta = expression(tau * sqrt(delta[t])),
          sigma = expression(sigma)
        )
      ) +
      labs(x = time, y = "SD induced by the time-varying covariate imbalance") +
      theme_bw()
    print(p)
  }
  delta
}
