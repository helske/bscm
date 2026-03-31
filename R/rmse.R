#' @export
#' @rdname rmse
rmse <- function(x, ...) {
  UseMethod("rmse", x)
}
#' Extract root mean squared errors of a Bayesian synthetic control model
#'
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param average \[`logical(1)`]\cr If `TRUE` (the default), returns the 
#' average RMSEs over treated units in case of multiple treated units. 
#' If `FALSE`, unit-specific values are returned.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the RMSE values.
#' @rdname rmse
#' @aliases rmse
#' @export
rmse.bscmfit <- function(x, probs = c(0.025, 0.975), average = TRUE, ...) {
  
  test_probs(probs)
  stopifnot_(
    checkmate::test_flag(average),
    "Argument {.arg average} must be a single {.cls logical} value."
  )
  N <- get_N(x)
  if (average && N > 1) {
    as_draws(x, c("avg_RMSE_pre", "avg_RMSE_post", "avg_RMSE_ratio")) |> 
      summarise_with_probs(probs) |>
      mutate(
        variable = case_when(
          startsWith(variable, "avg_RMSE_pre") ~ "Average pre-RMSE",
          startsWith(variable, "avg_RMSE_post") ~ "Average post-RMSE",
          startsWith(variable, "avg_RMSE_ratio") ~ "Average RMSE ratio"
        )
      )
  } else {
    unit <- get_unit(x)
    treated <- get_treated(x)
    as_draws(x, c("RMSE_pre", "RMSE_post", "RMSE_ratio")) |> 
      summarise_with_probs(probs) |>
      mutate("{unit}" := rep(treated, each = 3), .before = 1L) |> 
      mutate(
        variable = case_when(
          startsWith(variable, "RMSE_pre") ~ "Pre-RMSE",
          startsWith(variable, "RMSE_post") ~ "Post-RMSE",
          startsWith(variable, "RMSE_ratio") ~ "RMSE ratio"
        )
      )
  }
}
