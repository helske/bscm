
# bscm

<!-- badges: start -->
<!-- badges: end -->

R package for Bayesian model-based synthetic control method. Currently has only 
basic functionalities for inference in case of a single treated unit.

## Installation

You can install the development version of `bscm` as

``` r
remotes::install_github("helske/bscm")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bscm)
fit <- bscm(
  y ~ treatment, 
  data = simulated_example, 
  treatment = "treatment", 
  time = "time", 
  unit = "id"
)
fit$avg_effect # average treatment effects for pre- and post-treatment periods
fit$R2 # Bayesian R-squared estimate of the model fit
```

