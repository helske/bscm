
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
MCMC sampling time: 0.827 seconds

MCMC diagnostics indicate no issues. 
```
```{r}
> summary(fit)

  id    variable                   mean     sd   q2.5  q97.5  rhat ess_bulk ess_tail
1 1     Intercept                 2.46  0.571  1.27    3.51  1.00     4759.    3529.
2 NA    Coef_x                    1.05  0.429  0.214   1.91  1.000    4406.    2924.
3 1     Residual SD               0.521 0.0753 0.399   0.691 1.00     3748.    2729.
4 1     Bayesian R-squared        0.811 0.0484 0.700   0.886 1.00     4003.    2937.
5 NA    Average treatment effect  5.97  0.281  5.43    6.55  1.00     4124.    3851.
6 1     Pre-RMSE                  0.721 0.105  0.540   0.951 1.000    3993.    3770.
7 1     Post-RMSE                 6.69  0.306  6.11    7.32  1.00     4143.    3647.
8 1     RMSE ratio                0.108 0.0163 0.0794  0.143 1.000    3873.    3630.
9 1     Effective donors         11.7   2.35   7.10   16.3   1.00     2394.    3315.
```
And default visualization:
``` r
plot(fit)
```
![Figure showing synthetic control and treatment effect estimates](readme_fig.png)
