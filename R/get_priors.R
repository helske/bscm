#' @export
#' @rdname get_priors
get_priors <- function(x, ...) {
  UseMethod("get_priors")
}
#' Get prior specifications used in the Bayesian synthetic control model
#'
#' @param x \[`bscmfit`]\cr Output from [bscm()].
#' @param as_data_frame \[`logical(1)`]\cr If `TRUE` (the default),
#'   returns priors as a `data.frame`.
#' @param ... Ignored.
#' @return A named list of `bscm_prior` objects, or a `data.frame` with one
#'   row per element of each model parameter and columns `parameter`, 
#'   `term` (predictor), a column named by the unit variable (treated
#'   unit), and `prior`. Columns `term` and the unit column are omitted when
#'   they are not needed.
#' @export
#' @rdname get_priors
get_priors.bscmfit <- function(x, as_data_frame = TRUE, ...) {
  check_flag(as_data_frame, "as_data_frame")
  pr <- x$priors
  if (!as_data_frame) {
    return(pr)
  }
  weights <- pr$weights
  if (inherits(weights$concentration, "bscm_prior")) {
    pr$concentration <- weights$concentration
  }
  labels <- lapply(pr, as.character)
  if (!is.null(weights$tilt)) {
    labels$tilt <- paste0("normal(0, ", weights$tilt$scale, ")")
  }
  d <- lapply(
    names(labels),
    \(nm) {
      # the prior of the tilt coefficients is shared by units and directions
      elements <- if (nm == "tilt") {
        dplyr::tibble(term = NA, unit = NA)
      } else {
        parameter_elements(x, nm)
      }
      dplyr::tibble(parameter = nm, elements, prior = labels[[nm]])
    }
  ) |>
    dplyr::bind_rows() |>
    dplyr::rename("{get_unit(x)}" := "unit")
  as.data.frame(dplyr::select(d, dplyr::where(\(v) !all(is.na(v)))))
}
