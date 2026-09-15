# bscm 1.0.2

* The scale of the logistic normal prior of the donor weights `omega` can now
  be estimated from the data instead of being fixed by the user. Calling
  `logistic_normal_pr()` without arguments (equivalently,
  `logistic_normal_pr(NA)`) leaves the scale unspecified, in which case it is
  a unit-specific model parameter `zeta`. The prior for `zeta` is defined as
  the element `zeta` of the `priors` list of `bscm`.
* The elements of the `priors` list can now be named either by the name of the
  parameter in the posterior draws or by a plain-text synonym: `intercept`,
  `coef`, `scale_gamma`, `slope_gamma`, `ar1`, `weights`, `sd`, and `scale_omega` for `intercept`, `beta`, `kappa`, `nu`, `rho`, `omega`, `zeta`, and `sigma`
  respectively.
* The posterior summaries shown by `print` now use plain-text labels
  throughout, e.g. `scale_gamma_x` in place of `kappa_x` and `coef_x` in 
  place of `beta_x`. 
* Fixed `loo_R2` for multiple treated units.
* `coef` now returns the total effect of each predictor instead of the
  individual model parameters, one type at a time: the intercepts with
  `type = "intercept"`, the coefficients of the predictors whose effect is
  constant with `type = "fixed"` (the default), and the total effect
  `beta + gamma_t` of each time-varying predictor at each time point with
  `type = "varying"`.
* `plot_coefs` no longer has an argument `combine`, and always shows the total 
  effect `beta + gamma_t` of the time-varying predictors. 
* `summary` now returns posterior summaries of the model parameters, named as
  in the posterior draws, and no longer has a class and a print method of its
  own. It gained an argument `parameters` for selecting the parameters to
  summarise.
* `print` now shows the posterior summaries of the model parameters and the
  average pre- and post-treatment effects. With a single treated unit, the
  Bayesian R2 and the pre- and post-treatment RMSEs are shown as well, and
  with several treated units, only the parameters which are invariant over
  both time and unit.
* The simplex weight vector with Dirichlet prior is now based on the 
  isometric log-ratio transformation instead of stick-breaking transformation. 
  This leads considerably better behaving sampler with small concentration 
  parameters of the Dirichlet distribution. Note: The old and new constructions 
  correspond to pre-2.37 and post-2.37 versions of Stan, respectively.
* Due to improved stability of the new simplex construction, default 
 `adapt_delta` is now 0.8 for models without time-varying coefficients while the 
  old default 0.95 is still used with time-varying coefficients.
* Changes to the definition of the splines for time-varying 
  coefficients via `tv()` in order to make the pre-treatment coefficients 
  invariant to the length of the post-treatment period. For this purpose, 
  argument `df` now refers to the number of the basis functions in the 
  pre-treatment period using extended range of knot positions to avoid 
  boundary effects. In addition, `tv()` now has a new argument `knot_spacing` 
  which can be used as alternative to `df` to specify the number of knots in 
  the pre-treatment period.
* Default `type` argument in `tv()` is now `rw2` matching the second difference 
  penalty of the spline coefficients (old default was `rw1` for first 
  differences). Use of `rw2` adds another parameter `nu`, which is the slope 
  of the linear trend in `gamma`.
* Fixes and changes to `plot` methods: visualization is now separated from the 
  the summary methods, with argument `unit` specifying which units to plot.
* Fixed the generation of posterior predictive means for the post-treatment 
  period when using AR(1) residuals.
* Output of `loo` now contains an attribute `diagnostics` 
  with unit and time indices of observations with corresponding pareto-k 
  diagnostics.
* `proj_predict_bscm` has now arguments `summary` and `probs` for returning
  posterior predictive summaries. With `summary = FALSE`, the draws are now
  returned in the same long format as in e.g. `synthetic_control`, with the
  predictions as an `rvar` column. The function now also checks that the
  projection is based on a `bscmfit` object.
* Exported `get_treated` and `get_donors` method.
* Added `as_draws_rvars` method which returns the posterior draws as a
  `draws_rvars` object.
* Fixed `coef` with `summary = FALSE`, which silently dropped the
  autoregressive coefficients `rho` from the returned list.
* Fixed `covariate_adjustment`, which returned `NA` for predictors with 
  time-varying coefficients.
* Fixed `covariate_adjustment` and `covariate_imbalance` for factor and 
  interaction predictors, which previously failed with an error.
* Fixed the behaviour of `rmse_ratio` with `average = TRUE`: Now returns the 
  posterior of average ratio over units, instead of the ratio of averages.
* `residual_acf`, `covariate_imbalance` and `covariate_adjustment` now behave 
  like the other posterior summary methods with summary argument and any number 
  of quantiles in `probs`.
* The chain structure of the posterior draws is now retained when computing
  quantities that are based on the generated quantities `y_mean` and
  `y_rep`. Consequently the convergence diagnostics (`rhat`, `ess_bulk`,
  `ess_tail`, and `mcse_mean`) of `fitted`, `residuals`, `residual_acf`,
  `synthetic_control`, `treatment_effect`, `average_treatment_effect`,
  `rmse`, `rmse_ratio`, and `bayes_R2` are now computed across chains
  instead of treating all draws as a single chain. The posterior estimates
  themselves are unchanged.
* The `average` argument now defaults to `FALSE` in every method which has
  it. Previously `treatment_effect` and `average_treatment_effect` defaulted
  to `TRUE`, i.e. to averaging over the treated units.
* The unit column is no longer dropped from the output when the model
  contains a single treated unit. It is now absent only when the quantity was
  averaged over the treated units with `average = TRUE`, so that the shape of
  the output no longer depends on the number of treated units.
* Minor improvements to documentation and examples.

# bscm 1.0.1

* Initial CRAN submission.