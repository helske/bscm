
# bscm: Bayesian Synthetic Control Models

<!-- badges: start -->
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/helske/bscm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/helske/bscm/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/helske/bscm/graph/badge.svg)](https://app.codecov.io/gh/helske/bscm)
[![CRAN version](http://www.r-pkg.org/badges/version/bscm)](https://CRAN.R-project.org/package=bscm)
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
    treatment = "treatment",  time = "time",  unit = "id",
    chains = 2, refresh = 0, 
    priors = list(omega = dirichlet_pr(0.5), beta = normal_pr(0, 2))
)
```
Basic summary of the estimated model:
``` r
Call:
bscm(formula = y ~ x, data = single_treated, treatment = "treatment", 
    time = "time", unit = "id", priors = list(omega = dirichlet_pr(0.5), 
        beta = normal_pr(0, 2)), chains = 2, refresh = 0)

Bayesian synthetic control model y ~ x  
Treated unit: 1 
Number of donors: 50 
Number of time periods (pre + post): 40 + 10 
MCMC sampling using 2 chains, each with 2500 + 2500 iterations took 5.77 seconds for the slowest chain

MCMC diagnostics indicate no issues. 

Summary statistics of the time-invariant parameters:
  variable     mean     sd   q2.5 q97.5
1 Intercept   0.455 0.395  -0.313 1.23 
2 beta_x      1.01  0.0687  0.870 1.14 
3 Residual SD 0.678 0.0977  0.515 0.900

Average treatment effects and measures of model fit:
  variable                            mean      sd   q2.5 q97.5
1 Average pre-treatment effect  -0.0000702 0.156   -0.303 0.314
2 Average post-treatment effect  6.33      0.331    5.69  6.99 
3 Bayesian R2                    0.991     0.00271  0.984 0.995
4 Pre-treatment RMSE             0.946     0.136    0.714 1.24 
5 Post-treatment RMSE            7.29      0.340    6.64  7.99
```
And default visualization:
``` r
plot(fit)
```
![Figure showing synthetic control and treatment effect estimates](readme_fig.png)
