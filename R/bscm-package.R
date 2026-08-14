#' The 'bscm' package
#'
#' @description Implements the synthetic control method of Abadie, Diamond,
#' and Hainmueller (2010) within a Bayesian framework, enabling straightforward
#' uncertainty quantification of treatment effects and other quantities of
#' interest. Posterior sampling is performed using Markov chain Monte Carlo via
#' Stan.
#'
#' `bscm` supports time-varying covariates with potentially time-varying
#' effects, single treated unit as well as multiple treated units possibly with
#' staggered treatment adoption.
#'
#' Output of [bscm()] is a `bscmfit` object with methods for model assessment,
#' comparison, and selection based on placebo studies, cross-validation, and
#' posterior predictive checks.
#'
#' @name bscm-package
#' @useDynLib bscm, .registration = TRUE
#' @importFrom rstan sampling
#' @importFrom loo loo
#' @importFrom stats sigma coef fitted residuals sd
#' @importFrom projpred get_refmodel
#' @importFrom rstantools rstan_config log_lik bayes_R2 loo_R2
#' @importFrom rstantools posterior_epred posterior_linpred posterior_predict
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom methods new
#' @importFrom rlang := .data .env
#' @importFrom posterior as_draws nchains ndraws
#' @import ggplot2
#' @references
#' Abadie A, Diamond A, and Hainmueller J (2010). Synthetic Control
#' Methods for Comparative Case Studies: Estimating the Effect of California’s
#' Tobacco Control Program. *Journal of the American Statistical Association*,
#' 105(490), 493--505, <doi:10.1198/jasa.2009.ap08746>.
#'
#' Stan Development Team (2025). RStan: the R interface to Stan. R package
#' version 2.32.7. https://mc-stan.org.
#'
"_PACKAGE"
