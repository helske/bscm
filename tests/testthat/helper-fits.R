tryCatch(
  suppressWarnings({
    fit1_int <- bscm(
      y ~ 1,
      data = single_treated,
      treatment = "treatment",
      time = "time",
      unit = "id",
      chains = 2,
      iter = 2000,
      refresh = 0,
      seed = 1,
      save_warmup = FALSE
    )
    fit1_noint <- bscm(
      y ~ 0,
      data = single_treated,
      treatment = "treatment",
      time = "time",
      unit = "id",
      error = "ar1",
      chains = 2,
      iter = 200,
      init = 0,
      refresh = 0,
      seed = 1,
      save_warmup = FALSE
    )
    fit1_x <- bscm(
      y ~ x,
      data = single_treated,
      treatment = "treatment",
      time = "time",
      unit = "id",
      chains = 2,
      iter = 200,
      init = 0,
      refresh = 0,
      seed = 1,
      save_warmup = FALSE
    )
    fitN_int <- bscm(
      y ~ 1,
      data = multiple_treated,
      treatment = "treatment",
      time = "time",
      unit = "id",
      chains = 2,
      iter = 200,
      init = 0,
      refresh = 0,
      seed = 1,
      save_warmup = FALSE
    )
    fitN_xz <- bscm(
      y ~ x + z,
      data = multiple_treated,
      treatment = "treatment",
      time = "time",
      unit = "id",
      chains = 2,
      iter = 200,
      init = 0,
      refresh = 0,
      seed = 1,
      save_warmup = FALSE
    )
    fitN_tv <- bscm(
      y ~ z + tv(~x, 10),
      data = multiple_treated,
      treatment = "treatment",
      time = "time",
      unit = "id",
      priors = list(
        intercept = normal_pr(0, 1),
        kappa = exponential_pr(1)
      ),
      chains = 2,
      iter = 200,
      init = 0,
      refresh = 0,
      seed = 1,
      save_warmup = FALSE
    )
  }),
  error = function(e) {
    message("helper-fits.R: model fitting failed (", conditionMessage(e), ")")
  }
)
