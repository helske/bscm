#' Re-estimate bscmfit object x using pre-created standata
#'
#'@noRd
refit_bscm <- function(x, standata) {
  init <- lapply(x$stanfit@stan_args, `[[`, "init")
  stan_args <- x$stanfit@stan_args[[1]]
  stan_args <- list(
    object <- stanmodels$bscm,
    data <- standata,
    chains = length(init),
    iter = stan_args$iter,
    thin = stan_args$thin,
    warmup = stan_args$warmup,
    algorithm = stan_args$algorithm,
    control = stan_args$control,
    refresh = 0L,
    init = init
  )
  if (!is.null(x$setup$excluded_pars)) {
    stan_args$pars <- x$setup$excluded_pars
    stan_args$include <- FALSE
  }
  fit <- do.call(rstan::sampling, stan_args)

  gq <- rstan::gqs(
    stanmodels$generated_quantities,
    data = standata,
    draws = rstan::as.matrix(fit)
  )

  out <- list(
    stanfit = fit,
    y_mean = rstan::as.matrix(gq, "y_mean"),
    y_rep = "Not sampled",
    data = NULL,
    setup = x$setup,
    priors = x$priors,
    standata = standata
  )
  class(out) <- c("bscmfit_refit", "bscmfit")
  out
}
