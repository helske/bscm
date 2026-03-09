suppressWarnings({
  fit_intercept <- bscm(
    y ~ 1,
    data = simulated_data,
    treatment = "treatment",
    time = "time",
    unit = "id",
    chains = 1,
    iter = 300,
    warmup = 150,
    refresh = 0,
    seed = 123
  )
  fit_no_intercept <- bscm(
    y ~ 0,
    data = simulated_data,
    treatment = "treatment",
    time = "time",
    unit = "id",
    chains = 1,
    iter = 300,
    warmup = 150,
    refresh = 0,
    seed = 123
  )
  fit_x <- bscm(
    y ~ x,
    data = simulated_data,
    treatment = "treatment",
    time = "time",
    unit = "id",
    chains = 1,
    iter = 300,
    warmup = 150,
    refresh = 0,
    seed = 123
  )
  fit_tv_x <- bscm(
    y ~ x,
    time_varying_effects = TRUE,
    data = simulated_data,
    treatment = "treatment",
    time = "time",
    unit = "id",
    chains = 1,
    iter = 300,
    warmup = 150,
    refresh = 0,
    seed = 123,
  )
})
