# bscm 1.0.2

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
  differences).
* Fixes and changes to `plot` and `plot_effects` methods in case of multiple 
  treated units. Both methods now include argument `unit` to specify which 
  treated unit to plot. If `unit` is not specified, `plot_effects` plots 
  average treatment effect single treatment start, whereas `plot` plots all 
  units sequentially.
* Fixed the generation of posterior predictive means for the post-treatment 
  period when using AR(1) residuals.
* Output of `loo` now contains an attribute `diagnostics` 
  with unit and time indices of observations with corresponding pareto-k 
  diagnostics.
* `proj_predict_bscm` has now arguments `summary` and `probs` for returning
  posterior predictive summaries. 
* Exported `get_treated` and `get_donors` method.
* Minor improvements to documentation and examples.

# bscm 1.0.1

* Initial CRAN submission.