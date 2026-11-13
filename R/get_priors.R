#' @export
#' @rdname get_priors
get_priors <- function(x, ...) {
  UseMethod("get_priors")
}
#' Get prior specifications used in the Bayesian synthetic control model
#'
#' @param x \[`bscmfit`]\cr Output from [bscm()].
#' @param as_table \[`logical(1)`]\cr If `TRUE` (the default),
#'   returns priors as a `tibble`.
#' @param ... Ignored.
#' @return A named list of `bscm_prior` objects, or a `tibble` with one
#'   row per element of each model parameter and columns `parameter`, 
#'   `term` (predictor), a column named by the unit variable (treated
#'   unit), and `prior`. Column `term` is omitted when the model has no
#'   predictors, and `NA` values are printed as empty cells.
#' @export
#' @rdname get_priors
get_priors.bscmfit <- function(x, as_table = TRUE, ...) {
  check_flag(as_table, "as_table")
  pr <- x$priors
  if (!as_table) {
    return(pr)
  }
  weights <- pr$weights
  estimated <- inherits(weights$concentration, "bscm_prior")
  if (estimated) {
    pr$concentration <- weights$concentration
  }
  labels <- lapply(pr, as.character)
  # the estimated concentration and the tilt are shown on their own rows
  labels$weights <- if (estimated) {
    weights$distribution
  } else {
    paste0(weights$distribution, "(", weights$concentration, ")")
  }
  if (!is.null(weights$tilt)) {
    labels$tilt <- paste0(
      "normal(0, ", weights$tilt$scale, "), R = ", weights$tilt$R
    )
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
  if (all(is.na(d$term))) {
    d$term <- NULL
  }
  class(d) <- c("bscm_priors", class(d))
  d
}
#' Print the priors returned by `get_priors()` with empty cells for `NA`
#'
#' @param x \[`bscm_priors`]\cr Output of [get_priors()].
#' @param ... Further arguments passed to [print.data.frame()].
#' @noRd
#' @export
print.bscm_priors <- function(x, ...) {
  d <- x
  # a tibble would quote all values of a column containing empty strings
  class(d) <- "data.frame"
  d[is.na(d)] <- ""
  print(d, right = FALSE, ...)
  invisible(x)
}
