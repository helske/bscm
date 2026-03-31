
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
> fit
Call:
bscm(formula = y ~ x, data = single_treated, treatment = "treatment", 
    time = "time", unit = "id")

Bayesian synthetic control model y ~ x 
Treated unit: 1 
Number of donors: 30 
Number of time periods (pre + post): 30 + 10 
MCMC sampling time: 2.436 seconds

MCMC diagnostics indicate no major issues. 
```
```{r}
> summary(fit)

  id    variable                   mean     sd   q2.5  q97.5  rhat ess_bulk ess_tail
1 1     Intercept                 2.47  0.580  1.30    3.61  1.000    4824.    3743.
2 NA    Coef_x                    1.03  0.433  0.188   1.87  1.000    4683.    2986.
3 1     Residual SD               0.519 0.0741 0.397   0.683 1.00     3868.    2988.
4 1     Bayesian R-squared        0.813 0.0475 0.705   0.886 1.00     3895.    2469.
5 NA    Average treatment effect  5.96  0.281  5.42    6.52  1.000    4389.    3726.
6 1     Pre-RMSE                  0.720 0.102  0.543   0.946 1.00     3407.    3280.
7 1     Post-RMSE                 6.68  0.309  6.09    7.30  1.000    4142.    3760.
8 1     RMSE ratio                0.108 0.0160 0.0805  0.143 1.00     3458.    3360.
9 1     Effective donors         11.6   2.35   7.00   16.2   1.000    2656.    3496.
```
And default visualization:
``` r
plot(fit)
```
![Figure showing synthetic control and treatment effect estimates](readme_fig.png)
