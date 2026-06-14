#' Donor ranking
#'
#' Returns ranking of donors either using posterior means of donor weights from
#' `bscmfit` object, or ranking from a `vsel` object
#' returned by [projpred::varsel()] or [projpred::cv_varsel()] when applied to
#' a `bscmfit` object.
#'
#' @param x \[`vsel`] or \[`bscmfit`]\cr Output from [bscm()],
#' [projpred::varsel()] or [projpred::cv_varsel()].
#' @return Character vector of donor IDs in order of selection.
#' @export
donor_ranking <- function(x) {
  is_bscmfit <- inherits(x, "bscmfit")
  is_vsel <- inherits(x, "vsel")
  stopifnot_(
    is_bscmfit | is_vsel,
    "{.arg x} must be an object of class {.cls bscmfit} from {.fn bscm} or 
    class {.cls vsel} from {.fn projpred::varsel} or {.fn projpred::cv_varsel}."
  )
  if (is_bscmfit) {
    order_donors(x, order = "descending")
  } else {
    name_map <- x$refmodel$fit$proj$original_donor_names
    unname(name_map[projpred::ranking(x)$fulldata])
  }
}

#' Return donor ordering for leave one out and placebo runs
#' @noRd
order_donors <- function(x, order = NULL, weights = NULL) {
  donors <- get_donors(x)
  # original order
  if (is.null(order)) {
    return(donors)
  }
  # ranked ordering
  if (length(order) == 1L && order %in% c("ascending", "descending")) {
    if (is.null(weights)) {
      weights <- donor_weights(x, probs = numeric(0))
    }
    ranked <- weights |>
      dplyr::summarise(mean = mean(mean), .by = "donor") |>
      dplyr::arrange(if (order == "descending") dplyr::desc(mean) else mean)
    return(ranked$donor)
  }
  stopifnot_(
    checkmate::test_character(
      order,
      any.missing = FALSE,
      unique = TRUE,
      max.len = length(donors),
      min.len = 1L
    ),
    "Argument {.arg order} contains duplicate donor names."
  )
  stopifnot_(
    checkmate::test_subset(order, donors),
    "Argument {.arg order} contains names not corresponding to any donors."
  )
  order
}
