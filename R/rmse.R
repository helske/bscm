#' @export
#' @rdname rmse
rmse <- function(x, ...) {
  UseMethod("rmse", x)
}
#' Extract root mean squared errors of a Bayesian synthetic control model
#'
#' @export
#' @rdname rmse
#' @param x \[`bscmfit`]\cr The model fit object.
#' @param probs \[`numeric()`]\cr Probabilities for quantile summaries. 
#'   Default is `c(0.025, 0.975)`.
#' @param ... Ignored.
#' @return A `data.frame` of posterior summaries of the RMSE values.
rmse.bscmfit <- function(x, probs = c(0.025, 0.975), ...) {
  
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
  as_draws(x, c("RMSE_pre", "RMSE_post", "RMSE_ratio")) |> 
    summarize_with_probs(probs) |> 
    mutate(
      variable = case_when(
        variable == "RMSE_pre" ~ "Pre-RMSE",
        variable == "RMSE_post" ~ "Post-RMSE",
        variable == "RMSE_ratio" ~ "Post-RMSE / Pre-RMSE"
      )
    )
}
