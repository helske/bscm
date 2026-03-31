#' Visualize BSCM estimates
#' 
#' A plot of the posterior mean and posterior interval of the treatment effect 
#' and synthetic control over time.
#' 
#' @param x \[`bscmfit`]\cr object.
#' @param probs \[`numeric(2)`]\cr Vector of length two defining the limits of 
#' the posterior interval. Default is `c(0.025, 0.975)`.
#' @param ... Ignored
#' @aliases plot
#' @export
plot.bscmfit <- function(x, probs = c(0.025, 0.975), ...) {
  type <- yintercept <- ymin <- ymax <- NULL
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
  outcome <- get_outcome(x)
  treated <- get_treated(x)
  time <- get_time(x)
  unit <- get_unit(x)
  
  d_effects <- treatment_effect(x, probs, average = FALSE, for_plots = TRUE)
  d_synth <- synthetic_control(x, probs, for_plots = TRUE)
  d_effects$type <- "Treatment effect"
  d_synth$type <- "Synthetic control"
  lookup <- stats::setNames(
    c(unit, time, outcome, paste0("q", 100 * probs)), 
    c("unit", "time", "y", "ymin", "ymax"))
  d_plot <- bind_rows(d_effects, d_synth) |> 
    rename(any_of(lookup))
  
  dy <- x$data |> 
    filter(.data[[unit]] %in% .env$treated) |> 
    select(all_of(c(unit, time, outcome))) |> 
    mutate(type = "Synthetic control") |> 
    rename(any_of(c(mean = "y", lookup)))
  
  dt <- data.frame(yintercept = 0, type = "Treatment effect")
  N <- get_N(x)
  plots <- stats::setNames(vector("list", length = N), treated)
  for (i in treated) {
    plots[[i]] <- d_plot |> 
      filter(unit == i) |> 
      ggplot(aes(time, mean)) + 
      geom_hline(
        data = dt, aes(yintercept = yintercept), 
        linetype = "dashed", colour = "grey50"
      ) +
      geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = type), alpha = 0.25) +
      geom_line(aes(colour = type)) +
      geom_point(data = dy |> filter(unit == i), colour = "grey30") +
      labs(x = time, y = if (N > 1L) i) +
      scale_colour_manual(values = c("#DDAA33", "#0C7BDC")) +
      scale_fill_manual(values = c("#EECC66", "#77AADD")) +
      guides(fill = "none", colour = "none") +
      facet_wrap(~ type, ncol = 1, scales = "free_y", strip.position = "left") +
      theme_bw()
  }
  if (N == 1L) plots[[1]] else plots
}