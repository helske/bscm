
# bscm: Bayesian Synthetic Control Models

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/helske/bscm/graph/badge.svg)](https://app.codecov.io/gh/helske/bscm)
[![R-CMD-check](https://github.com/helske/bscm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/helske/bscm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

R package for Bayesian Synthetic Control Models (Helske 2026, in-preparation). 
Key features:
* Time-varying covariates, optionally with time-varying coefficients as splines.
* Multiple treated units, with or without staggered adoption.
* Weakly informative default priors and user-defined priors are supported.
* Computationally and statistically efficient posterior sampling via 
  pre-compiled [Stan](https://mc-stan.org/) models.
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
MCMC sampling using 4 chains, each with 2500 + 2500 iterations took 6.82 seconds for the slowest chain

MCMC diagnostics indicate no issues. 
```
``` r
summary(fit)
  variable                          mean      sd   q2.5  q97.5  rhat ess_bulk ess_tail mcse_mean
1 Intercept                      0.506   0.370   -0.223  1.24   1.00   11183.    9289. 0.00350  
2 beta_x                         1.01    0.0703   0.874  1.15   1.00   10267.    7509. 0.000695 
3 Residual SD                    0.781   0.112    0.595  1.03   1.00    8848.    8223. 0.00120  
4 Bayesian R2                    0.988   0.00352  0.979  0.993  1.00    8757.    8100. 0.0000384
5 Effective number of donors    23.0     3.70    15.2   29.8    1.00    6816.    8774. 0.0448   
6 Average pre-treatment effect  -0.00158 0.175   -0.352  0.339  1.00   10068.    8911. 0.00175  
7 Average post-treatment effect  6.57    0.331    5.94   7.24   1.00   10120.    9078. 0.00329  
8 Pre-treatment RMSE             1.08    0.152    0.821  1.41   1.00    8384.    8599. 0.00167  
9 Post-treatment RMSE            7.61    0.332    6.99   8.29   1.00    9860.    9234. 0.00335  
```
And default visualization:
``` r
plot(fit)
```
![Figure showing synthetic control and treatment effect estimates](readme_fig.png)
