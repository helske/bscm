#' @export
#' @rdname get_priors
get_priors <- function(x, ...) {
  UseMethod("get_priors")
}
#' Get prior specifications used in the Bayesian synthetic control model
#' @param x \[`bscmfit`]\cr Output from [bscm()].
#' @param as_data_frame \[`logical(1)`]\cr If `TRUE`,
#' returns priors as a `data.frame`. Default is `TRUE`.
#' @param ... Ignored.
#' @return A named list of `bscm_prior` objects for the parameters of the model.
#' @export
#' @rdname get_priors
get_priors.bscmfit <- function(x, as_data_frame = TRUE, ...) {
  check_flag(as_data_frame, "as_data_frame")
  pr <- x$priors
  if (as_data_frame) {
    pr <- data.frame(
      parameter = names(pr),
      distribution = unlist(lapply(pr, as.character), use.names = FALSE)
    )
  }
  pr
}
