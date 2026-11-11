# bscm 1.1.0

* The model parameters are now named by their roles in all user-facing
  input and output: `intercept`, `coef`, `tv_sd`, `tv_slope`, `sigma`, `ar1`,
  `weights`, `concentration`, and `tilt`, see `?bscm_prior` for their
  definitions and the corresponding names of the posterior draws (`alpha`,
  `beta`, `kappa`, `nu`, `sigma`, `rho`, `omega`, `zeta`, and `phi`).
* The standard deviations `sigma` of the error term are now integrated out
  analytically during MCMC sampling by default, and their posterior draws are
  sampled afterwards from their conditional posterior distribution. This is
  controlled by the new argument `marginalize_sigma` of `bscm`, and is
  supported with the priors `half_normal_pr()`, `scaled_inv_chi_pr()`, and
  `sqrt_inv_gamma_pr()` of `sigma`. The default prior of `sigma` is now
  `half_normal_pr(sd_e)`.
* Logistic normal prior (and new logistic Laplace prior) of donor weights is 
  now parameterized by concentration `zeta = 1/scale^2`, where scale is the 
  standard deviation of the log-weights.
* Concentration parameter of donor weights can now
  be estimated instead of being fixed, by defining argument `concentration` of 
  `logistic_normal_pr()`, `logistic_laplace_pr()`, and `dirichlet_pr()`, as 
  `bscm_prior, e.g., `logistic_normal_pr(inv_gamma_pr(10, 2))`. Calling these 
  functions without `concentration` uses the default prior, while a 
  number fixes it as before.
* The priors of the donor weights gained an argument `tilt`, which
  adds unit-specific coefficients `phi` that tilt the prior along the first 
  `R` principal components of the donors' pre-treatment outcomes 
  with `normal_pr(0, scale)` prior for `phi`.
* The default prior of the donor weights is now the logistic normal
  prior with estimated concentration `zeta` and a tilt along the first three
  principal components of the donors,
  `logistic_normal_pr(tilt = list(R = 3, scale = 2))`, instead of the uniform
  Dirichlet prior. Without marginalization of `sigma`, the default prior has
  no tilt.
* Added more prior distributions.
* `get_priors` now returns a data frame with columns `parameter`, `term`,
  the unit variable, and `prior`, with one row per predictor or treated unit,
  and the priors of the concentrations and the tilt coefficients on their own
  rows.
* Fixed `loo_R2` for multiple treated units.
* `coef` now returns the total effect of each predictor instead of the
  individual model parameters, one type at a time: the intercepts with
  `type = "intercept"`, the coefficients of the predictors whose effect is
  constant with `type = "fixed"` (the default), and the total effect
  `beta + gamma_t` of each time-varying predictor at each time point with
  `type = "varying"`. The output has columns `parameter`, `term` (the
  predictor), and the unit or time variable.
* `plot_coefs` no longer has an argument `combine`, and always shows the total 
  effect `beta + gamma_t` of the time-varying predictors. 
* `summary` now returns posterior summaries of the model parameters in long
  format, with columns `parameter` (role name), `term` (predictor or tilt
  direction), and the unit variable, and no longer has a class and a print
  method of its own. It gained an argument `parameters` for selecting the
  parameters to summarise by their role names.
* The column `treated` of the output of `donor_weights` is now named by the
  unit variable, and the column `variable` of `covariate_adjustment` is now
  `term`.
* All methods returning posterior summaries, such as `coef` and `rmse`, gain 
  new argument `diagnostics` to control whether to compute and return the MCMC 
  diagnostics such as Rhat and ESS values.
* `print` now shows the posterior summaries of the model parameters in
  sections "Coefficients", "Time variation of coefficients", "Error", and
  "Donor weights", with short labels such as `x`, `sd(x)`, and `slope(x)`,
  and the average pre- and post-treatment effects. With a single treated unit, 
  the Bayesian R2 and the pre- and post-treatment RMSEs are shown as well, and
  with several treated units, only the parameters which are invariant over
  both time and unit.
* The simplex weight vector with Dirichlet prior is now based on the 
  isometric log-ratio transformation instead of stick-breaking transformation. 
  This leads considerably better behaving sampler with small concentration 
  parameters of the Dirichlet distribution. Note: The old and new constructions 
  correspond to pre-2.37 and post-2.37 versions of Stan, respectively.
* Changes to the definition of the splines for time-varying 
  coefficients via `tv()` in order to make the pre-treatment coefficients 
  invariant to the length of the post-treatment period. For this purpose, 
  argument `df` now refers to the number of the basis functions in the 
  pre-treatment period using extended range of knot positions to avoid 
  boundary effects. In addition, `tv()` now has a new argument `knot_spacing` 
  as as alternative to `df` to specify the number of knots.
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