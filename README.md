
# bscm

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/helske/bscm/graph/badge.svg)](https://app.codecov.io/gh/helske/bscm)
[![R-CMD-check](https://github.com/helske/bscm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/helske/bscm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

R package for Bayesian Synthetic Control Models (Helske 2026, in-preparation). 


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
MCMC sampling using 4 chains, each with 2500 + 2500 iterations took 1.54 seconds

MCMC diagnostics indicate no issues. 
```
```{r}
summary(fit)
  id    variable                   mean     sd    q2.5  q97.5  rhat ess_bulk ess_tail
1 NA    Coef_x                0.995     0.447   0.0884  1.84  1.000    9158.    7341.
2 1     Residual SD (sigma)   0.512     0.0781  0.384   0.690 1.00     9895.    6800.
3 1     Bayesian R-squared    0.818     0.0499  0.699   0.894 1.00     9740.    7562.
4 1     Pre-treatment effect  0.0000376 0.133  -0.269   0.259 1.000   11252.    8803.
5 1     Post-treatment effect 5.92      0.300   5.34    6.51  1.000    9545.    9055.
6 1     Pre-RMSE              0.711     0.106   0.528   0.936 1.00     8220.    8255.
7 1     Post-RMSE             6.64      0.322   6.03    7.29  1.000    9248.    8696.
8 1     RMSE ratio            9.54      1.45    6.97   12.7   1.00     8437.    8577.
9 1     Effective donors      8.40      3.13    2.81   14.7   1.00     8730.    8809.
```
And default visualization:
``` r
plot(fit)
```
![Figure showing synthetic control and treatment effect estimates](readme_fig.png)
