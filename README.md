
# bscm: Bayesian Synthetic Control Models

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/helske/bscm/graph/badge.svg)](https://app.codecov.io/gh/helske/bscm)
[![R-CMD-check](https://github.com/helske/bscm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/helske/bscm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

R package for Bayesian Synthetic Control Models (Helske 2026, in-preparation). 
Key features:
* Time-varying covariates, optionally with time-varying random-walk coefficients.
* Multiple treated units, with or without staggered adoption.
* Customizable Dirichlet and logistic normal priors for donor weights.
* Convenient methods for extracting posterior summaries or draws of 
  treatment effects, synthetic control series, donor weights, RMSE, 
  Bayesian R^2, and other quantities of interest.
* Model evaluation and comparison via leave-one-out and leave-future-out 
  cross-validation, leave-donor(s)-out, in-time, and in-space placebo studies.


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
Number of donors: 30 
Number of time periods (pre + post): 30 + 10 
MCMC sampling using 4 chains, each with 2500 + 2500 iterations took 1.34 seconds for the slowest chain

MCMC diagnostics indicate no issues.  
```
```{r}
summary(fit)
  id    variable                    mean     sd   q2.5  q97.5  rhat ess_bulk ess_tail
1 NA    Coef_x                 1.04      0.413   0.223  1.86  1.00     9157.    7162.
2 1     Residual SD (sigma)    0.522     0.0748  0.399  0.692 1.000    9185.    7201.
3 1     Bayesian R-squared     0.809     0.0473  0.698  0.884 1.000    9797.    7824.
4 1     Pre-treatment effect  -0.0000958 0.137  -0.271  0.274 1.00    10145.    8758.
5 1     Post-treatment effect  5.97      0.261   5.45   6.49  1.00    10074.    9368.
6 1     Pre-RMSE               0.724     0.103   0.546  0.947 1.000    8832.    8503.
7 1     Post-RMSE              6.69      0.285   6.13   7.26  1.000    9421.    9228.
8 1     RMSE ratio             9.43      1.38    6.95  12.4   1.00     8936.    8927.
9 1     Effective donors      16.2       2.30   11.5   20.5   1.00     5244.    7033.
```
And default visualization:
``` r
plot(fit)
```
![Figure showing synthetic control and treatment effect estimates](readme_fig.png)
