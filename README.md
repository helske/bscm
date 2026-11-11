
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
    chains = 2, refresh = 0
)
```
Basic summary of the estimated model:
``` r
Call:
bscm(formula = y ~ x, data = single_treated, treatment = "treatment", 
    time = "time", unit = "id", chains = 2, refresh = 0)

Bayesian synthetic control model y ~ x  
Treated unit: 1 
Number of donors: 50 
Number of time periods (pre + post): 40 + 10 
MCMC sampling using 2 chains with 2500 + 2500 iterations took 12.5 seconds

MCMC diagnostics indicate no issues. 

Coefficients:
  variable   mean     sd   q2.5 q97.5
1 intercept 0.468 0.459  -0.421  1.38
2 x         0.936 0.0860  0.770  1.11

Error:
  variable  mean     sd  q2.5 q97.5
1 sigma    0.594 0.0766 0.464 0.761

Donor weights:
  variable        mean    sd   q2.5 q97.5
1 concentration  0.277 0.129  0.123 0.614
2 tilt(PC1)     -0.435 0.975 -2.53  1.37 
3 tilt(PC2)     -0.973 0.942 -3.01  0.734
4 tilt(PC3)      2.27  0.964  0.471 4.38 

Average treatment effects and measures of model fit:
  variable                          mean      sd   q2.5 q97.5
1 Average pre-treatment effect  0.000539 0.133   -0.266 0.266
2 Average post-treatment effect 5.50     0.443    4.67  6.41 
3 Bayesian R2                   0.993    0.00185  0.989 0.996
4 Pre-treatment RMSE            0.823    0.104    0.644 1.05 
5 Post-treatment RMSE           6.28     0.453    5.42  7.19 
```
And default visualization:
``` r
plot(fit)
```
![Figure showing synthetic control and treatment effect estimates](readme_fig.png)
