#' The 'bscm' package
#'
#' @description Bayesian synthetic control models for causal inference. Based on a
#' Bayesian model-based version of classic synthetic control approach.
#' Supports single or multiple treated units with staggered treatment times
#' and automatic uncertainty quantification based on posterior samples. The
#' fast estimation uses Markov chain Monte Carlo via Stan.
#'
#' @name bscm-package
#' @useDynLib bscm, .registration = TRUE
#' @importFrom rstan sampling extract
#' @importFrom loo loo
#' @importFrom stats sigma coef
#' @importFrom projpred get_refmodel
#' @importFrom rstantools rstan_config log_lik bayes_R2 loo_R2
#' @importFrom rstantools posterior_epred posterior_linpred posterior_predict
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom methods new
#' @importFrom rlang := .data .env
#' @import ggplot2
#' @import posterior
#' @import dplyr
#' @references
#' Stan Development Team (2025). RStan: the R interface to Stan. R package version 2.32.7. https://mc-stan.org
#'
"_PACKAGE"
