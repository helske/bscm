# bscm 1.0.2

* Fixes and changes to `plot` and `plot_effects` methods in case of multiple 
  treated units. Both methods now include argument `unit` to specify which 
  treated unit to plot. If `unit` is not specified, `plot_effects` plots 
  average treatment effect single treatment start, whereas `plot` plots all 
  units sequentially.
* Fixed the generation of posterior predictive means for the post-treatment 
  period when using AR(1) residuals.
* Output of `loo` with `reloo = TRUE` now contains an attribute `diagnostics` 
  with unit and time indices of observations with corresponding pareto-k 
  diagnostics.

# bscm 1.0.1

* Initial CRAN submission.