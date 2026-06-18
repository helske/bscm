
# bscm: Bayesian Synthetic Control Models

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/helske/bscm/graph/badge.svg)](https://app.codecov.io/gh/helske/bscm)
[![R-CMD-check](https://github.com/helske/bscm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/helske/bscm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

R package for Bayesian Synthetic Control Models (Helske 2026, in-preparation). 
Key features:
* Time-varying covariates, optionally with time-varying coefficients as splines.
* Multiple treated units, with or without staggered adoption.
* Customizable Dirichlet and logistic normal priors for donor weights.
* Convenient methods for extracting posterior summaries or draws of 
  treatment effects, synthetic control series, donor weights, RMSE, 
  Bayesian R^2, and other quantities of interest.
* Model evaluation and comparison using leave-one-out and leave-future-out 
  cross-validation, leave-donor(s)-out, in-time and in-space placebos.


## Installation

You can install the development version of `bscm` as

``` r
remotes::install_github("helske/bscm")
```

## Example

``` r
library(bscm)
set.seed(3546)
fit <- bscm(
    y ~ x, data = single_treated, 
    treatment = "treatment",  time = "time",  unit = "id"
)
```
Basic summary of the estimated model:
``` r
fit
Call:
bscm(formula = y ~ x, data = single_treated, treatment = "treatment", 
    time = "time", unit = "id")

Bayesian synthetic control model y ~ x  
Treated unit: 1 
Number of donors: 50 
Number of time periods (pre + post): 40 + 10 
MCMC sampling using 4 chains, each with 2500 + 2500 iterations took 3.9 seconds for the slowest chain

MCMC diagnostics indicate no issues. 
```
``` r
summary(fit)
  variable                  mean     sd   q2.5 q97.5  rhat ess_bulk ess_tail mcse_mean
1 Intercept              0.501   0.369  -0.220 1.22  1.00    12663.    9499.  0.00328 
2 beta_x                 1.01    0.0709  0.873 1.15  1.00    13039.    8017.  0.000620
3 Residual SD            0.782   0.109   0.600 1.02  1.00    12359.    8052.  0.000997
4 Pre-treatment effect  -0.00153 0.176  -0.345 0.346 1.000   12157.    8905.  0.00160 
5 Post-treatment effect  6.57    0.333   5.93  7.24  1.00    11150.    8971.  0.00316 
6 Pre-treatment RMSE     1.09    0.152   0.823 1.41  1.00    10114.    9050.  0.00152 
7 Post-treatment RMSE    7.62    0.336   6.97  8.30  1.00    11346.    9363.  0.00316 
```
And default visualization:
``` r
plot(fit)
```
![Figure showing synthetic control and treatment effect estimates](readme_fig.png)
