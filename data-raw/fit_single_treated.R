## code to used to create `fit_single_treated` dataset
set.seed(346)
fit_single_treated <- bscm(
  y ~ x,
  data = single_treated,
  treatment = "treatment",
  chains = 2,
  cores = 1,
  refresh = 0,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.8)
)
usethis::use_data(fit_single_treated, overwrite = TRUE)
