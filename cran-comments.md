## Comments for resubmission

* There are (not yet) any references for the proposed method. 
  These will be added in later updates.
* Added missing return values of functions in docs.
* Replaced `\dontrun{}` with `\donttest{}`.
* On modifying the global environment: 
  `loo_R2.bscmfit` function in `R/R2.R` had lines
  ```
  old_seed <- .Random.seed
  on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
  set.seed(get_stanfit(object)@stan_args[[1]]$seed)
  ```
  This was specifically done to avoid altering the global `Random.seed` by 
  setting and restoring its state on exit. This is now wrapped in if clause
  `if (fixed_seed)` where `fixed_seed` is argument to `loo_R2` 
  (with default `TRUE`).