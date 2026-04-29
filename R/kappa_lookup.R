#' Lookup table of weight vector priors and corresponding effective number of donors
#'
#' Based on Monte Carlo simulations, this data contains approximate 5%, 50%, and 
#' 95% quantiles of effective number of donors \eqn{1/\sum(\omega^2)} statistic 
#' for different values of \eqn{\kappa} parameter of Dirichlet and logistic normal 
#' priors of donor weight vectors of varying lengths. To see exactly how the data 
#' was generated, see `data-raw` folder on the Github repository of the package.
#'
#' @docType data
#' @keywords datasets
#' @format
#' A data frame with 4,350 rows and 9 columns:
#'  * distribution: Prior distribution.
#'  * J: Number of donors (length of weight vector)
#'  * kappa: Depending on the distribution, either concentration parameter 
#'    of symmetric Dirichlet distribution, or the scale parameter of 
#'    logistic normal distribution.
#'  * median_ess, q5_ess, q95_ess: Median, 5% and 95% quantiles of ESS.
#'  * median_r_ess, q5_r_ess, q95_r_ess: Median, 5% and 95% quantiles of 
#'   relative ESS i.e. ESS / J.
#' @name kappa_lookup
#' @examples
#' head(kappa_lookup)
NULL
