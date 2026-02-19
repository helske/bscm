#' Visualize BSCM estimates
#' 
#' A plot of the posterior mean and posterior interval of the treatment effect 
#' and synthetic control over time.
#' 
#' @param x \[`bscmfit`]\cr object.
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of 
#' the posterior interval. Default is `c(0.025, 0.975)`.
#' @param ... Ignored
#' @export
plot.bscmfit <- function(x, probs = c(0.025, 0.975), ...) {
  type <- yintercept <- ymin <- ymax <- NULL
  stopifnot_(
    checkmate::test_numeric(
      x = probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      len = 2L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector of length 2 with 
    values between 0 and 1."
  )
  outcome <- get_outcome(x)
  treated <- get_treated(x)
  time <- get_time(x)
  unit <- get_unit(x)
  T_pre <- get_T_pre(x)
  times <- get_times(x)
  
  draws <- as_draws(x, c("effect", "synthetic_y"))
  
  d_plot <- dplyr::bind_rows(
    "Synthetic control" = draws |> 
      subset_draws("synthetic_y") |> 
      summarise_draws(mean, ~ quantile2(.x, probs = probs)),
    "Treatment effect" = draws |> 
      subset_draws("effect") |> 
      summarise_draws(mean, ~ quantile2(.x, probs = probs)),
    .id = "type"
  ) |>
    mutate(time = rep(.env$times, 2), .before = 2L) |> 
    select(-"variable")
  colnames(d_plot) <- c("type", "time", "mean", "ymin", "ymax")
  
  y <- x$data |> 
    filter(.data[[unit]] %in% .env$treated) |> 
    pull(.data[[outcome]])
  
  dt <- data.frame(yintercept = 0, type = "Treatment effect")
  dy <- data.frame(time = times, mean = y, type = "Synthetic control")
  d_plot |> ggplot(aes(time, mean, group = type)) + 
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(
      data = dt, aes(yintercept = yintercept), 
      linetype = "dashed", colour = "grey50"
    ) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = type), alpha = 0.5) +
    geom_line(aes(colour = type), linewidth = 1) +
    geom_point(data = dy, colour = "grey30") +
    labs(x = time, y = NULL) +
    scale_colour_manual(values = c("#DDAA33", "#0C7BDC")) +
    scale_fill_manual(values = c("#EECC66", "#77AADD")) +
    guides(fill = "none", colour = "none") +
    facet_wrap(~ type, scales = "free_y", strip.position = "left", ncol = 1) +
    theme_bw()
}

