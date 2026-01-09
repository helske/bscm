#' Estimate Bayesian synthetic control model
#' 
#' Function `bscm` estimates a Bayesian synthetic control models of varying types.
#' 
#' Priors for intercept \eqn{\alpha} and standard deviation \eqn{\sigma} of the 
#' noise term are 
#' \eqn{\alpha \sim N(m, s)} and \eqn{\sigma \sim exponential(l)},
#' where by default `m = mean(y)`, `s = 2 * sd(y)`, `l = 1 / sd(diff(y))`,
#' with `y` being the vector of pretreatment outcomes. 
#' For multiple treated units, each unit have separate \eqn{\alpha} and 
#' \eqn{\sigma} and their prior definitions.
#' 
#' The prior for the weight vector \eqn{w} is defined as a symmetric Dirichlet 
#' distribution with concentration parameter \eqn{\kappa}. By default, 
#' \eqn{\kappa = 1}, corresponding to uniform distribution on a simplex. 
#' Alternatively, argument `effective_donors` can be used to define 
#' prior for \eqn{\kappa} using the approximation of expected number of 
#' effective donors, defined as \eqn{1 / sum(w^2)}. More specifically, 
#' `kappa = (effective_donors - 1) / (J - effective_donors)` where J is 
#' the total number of donors.
#' 
#' @param formula \[`formula`]\cr The model formula containing the outcome 
#' variable on the left-hand side and treatment variable and optional additional 
#' time-varying terms on the right-hand side of `~`.
#' @param data \[`data.frame` or an object coercible to one]\cr
#'   The long format data that contains the model variables.
#' @param treatment  \[`character(1)`]\cr Name of the treatment indicator 
#' variable in `data`.
#' @param time \[`character(1)`]\cr Name of the time index variable in `data`.
#' @param unit \[`character(1)`]\cr Name of the variable in `data` identifying 
#' different units.
#' @param priors \[`list()` or `character(1)`]\cr List of prior definitions 
#' or `"default"` which is a default and only supported option at the moment. 
#' See details on the prior definitions.
#' @param priors Definition of priors for intercept and standard deviation 
#' parameters. Currently, only the option `"default"` is supported. See details.
#' @param intercept Should the model include intercept term? Default is `TRUE`.
#' @param effective_donors Optional argument for defining a prior for weights. 
#' See details.
#' @param probs Vector of probabilities used for computing posterior quantiles.
#'  Defaults to `c(0.05, 0.95)`, i.e., the results contain limits of 90% 
#'  posterior intervals. Note that `probs` can contain more any number of values,
#'  e.g., `probs = c(0.025, 0.05, 0.5, 0.7, 0.9)` is valid.
#' @param ... Additional parameters passed on to [rstan::sampling()], e.g.,
#' `cores` and `iter`.
#' @return A list of data frames containing the estimates, as well as the full 
#' `stanfit` object, see [rstan::sampling()] for details.
#' @export
#' @examples
#' fit <- bscm(
#'   y ~ treatment, simulated_example, "treatment", "time", "id", 
#'   effective_donors = 10, chains = 1, cores = 1, refresh = 0
#' )
#' fit$avg_effect
#' fit$R2
#' 
bscm <- function(formula, data, treatment, time = "time", unit = "id",
                 priors = "default", intercept = TRUE,
                 effective_donors = NULL, probs = c(0.05, 0.95), ...) {
  check_bscm_arguments(
    formula, data, treatment, time, unit, priors, intercept, effective_donors, 
    probs
  )
  
  data <- data |>
    arrange(.data[[unit]], .data[[time]])
  outcome <- get_outcome(formula)
  treatment_table <- table(data[[treatment]], data[[unit]])
  treated <- names(which(treatment_table[2, ] > 0))
  controls <- names(which(treatment_table[2, ] == 0))
  T_pre <- treatment_table[1, treated] # TODO, change for multiple treated units
  T_total <- length(unique(data[[time]]))
  Y <- data |> 
    filter(.data[[unit]] %in% .env$treated) |> 
    pull(.data[[outcome]]) |> 
    matrix(nrow = T_total)
  Z <- data |> 
    filter(.data[[unit]] %in% .env$controls) |> 
    pull(.data[[outcome]]) |> 
    matrix(nrow = T_total)
  stopifnot_(
    all(!is.na(Y)) && all(!is.na(Z)),
    "Missing values are not (yet) supported."
  )
  J <- ncol(Z)
  f <- update.formula(formula, as.formula(paste0("~ . -", treatment)))
  predictors <- !identical(deparse(update(f, 0 ~ .)), "0 ~ 1")
  
  N <- length(treated)
  stopifnot_(
    N == 1,
    "Only the case of a single treated unit is currently implemented."
  )
  
  if (is.null(effective_donors)) {
    kappa <- 1 # uniform prior for weights
  } else {
    stopifnot_(
      effective_donors >= 2 && effective_donors < min(T_pre, J - 1),
      "Argument {.arg effective_donors} should be between 2 and 
      {min(T_pre, J - 1)}."
    )
    kappa <- (effective_donors - 1) / (J - effective_donors)
  }
  model_type <- case_when(
    N == 1 & !predictors & intercept ~ "single_nox",
    N == 1 & !predictors & !intercept ~ "single_nox_noalpha",
    N == 1 & predictors ~ "single_x",
    N > 1 & !predictors ~ "multiple_nox",
    N > 1 & predictors ~ "multiple_x",
  )
  single_nox_standata <- function(T_total, T_pre, Y, Z) {
    mean_y <- mean(Y[1:T_pre])
    sd_y <- sd(diff(Y[1:T_pre]))
    stopifnot_(
      sd_y > 0,
      "Outcome variable cannot be constant in the pretreatment period."
    )
    list(
      T = T_total, T_pre = T_pre, J = J, y = c(Y), Z = Z,
      pr_scale_sigma = 1 / sd_y, pr_mean_intercept = mean_y,
      pr_sd_intercept = 2 * sd_y, kappa = kappa
    )
  }
  standata <- case_when(
    model_type == "single_nox" ~ single_nox_standata(T_total, T_pre, Y, Z),
    model_type == "single_nox_noalpha" ~ single_nox_standata(T_total, T_pre, Y, Z),
    .default = NA
  )
  stopifnot_(
    is.list(standata),
    "Only single treated, no covariate case is currently implemented."
  )
  stan_args <- list(...)
  if (is.null(stan_args$init)) {
    if (is.null(stan_args$chains)) stan_args$chains <- 4L
    stan_args$init <- replicate(
      stan_args$chains, 
      list(
        omega = rep(1 / ncol(Z), ncol(Z)),
        alpha = rnorm(1, mean(Y[1:T_pre]), 0.1),
        sigma = runif(1, 0.5 * sd(Y[1:T_pre]), 2 * sd(Y[1:T_pre]))
      ),
      simplify = FALSE
    )
  }
  stan_args$object <- stanmodels[[model_type]]
  stan_args$data <- standata
  fit <- do.call(sampling, stan_args)
  draws <- as_draws(fit)
  times <- unique(data[[time]])
  if (N > 1) {
    stop("Multiple treated units are not yet supported.")
  } else {
    effect <- draws |> subset_draws("effect") |>
      summarise_draws(
        mean, sd, ~ quantile2(.x, probs = probs), default_convergence_measures()
      ) |> 
      mutate(time = times) |> 
      select(-variable) |> 
      relocate(time)
    colnames(effect)[1] <- time
    
    synthetic <- draws |> subset_draws("synthetic_y") |> 
      summarise_draws(
        mean, sd, ~ quantile2(.x, probs = probs), default_convergence_measures()
      ) |> 
      mutate(time = times) |> 
      select(-variable) |> 
      relocate(time)
    colnames(synthetic)[1] <- time
    
    weights <- draws |> subset_draws("omega") |> 
      summarise_draws(
        mean, sd, ~ quantile2(.x, probs = probs), default_convergence_measures()
      ) |> 
      mutate(donor = controls) |> 
      select(-variable) |> 
      relocate(donor)
    colnames(weights)[1] <- unit
    
    avg_eff <- draws |> subset_draws(c("avg_effect_pre", "avg_effect_post")) |> 
      summarise_draws(
        mean, sd, ~ quantile2(.x, probs = probs), default_convergence_measures()
      ) |> 
      mutate(period = factor(c("Pre-treatment", "Post-treatment"))) |> 
      select(-variable) |> 
      relocate(period)
    
    avg_cumeff <- draws |> subset_draws("avg_effect_post_cumulative") |> 
      summarise_draws(
        mean, sd, ~ quantile2(.x, probs = probs), default_convergence_measures()
      ) |> 
      mutate(time = times[(T_pre + 1):T_total]) |> 
      select(-variable) |> 
      relocate(time)
    colnames(avg_cumeff)[1] <- time
    
    rel_change <- draws |> subset_draws("relative_change") |> 
      summarise_draws(
        mean, sd, ~ quantile2(.x, probs = probs), default_convergence_measures()
      ) |> 
      mutate(time = times[(T_pre + 1):T_total]) |> 
      select(-variable) |> 
      relocate(time)
    colnames(rel_change)[1] <- time
    rmse <- draws |> subset_draws(c("pre_RMSE", "post_RMSE")) |> 
      summarise_draws(
        mean, sd, ~ quantile2(.x, probs = probs), default_convergence_measures()
      ) |> 
      mutate(period = factor(c("Pre-treatment", "Post-treatment"))) |> 
      select(-variable) |> 
      relocate(period)
    
    # RMSE of point estimates (posterior mean)
    rmse_pe <- effect |> 
      mutate(period = ifelse(row_number() <= T_pre, "Pre", "Post")) |> 
      summarise(RMSE = sqrt(mean(mean^2)), .by = period) |> 
      arrange(desc(period)) |> 
      pull(RMSE)
    rmse$RMSE_point <- rmse_pe
    
    R2 <- draws |> subset_draws("R2") |> 
      summarise_draws(
        mean, sd, ~ quantile2(.x, probs = probs), default_convergence_measures()
      ) |> 
      select(-variable)
    
    out <- list(
      effect = effect, synthetic = synthetic, weights = weights, 
      avg_effect = avg_eff, avg_cumulative_eff = avg_cumeff, 
      relative_change = rel_change, RMSE = rmse, R2 = R2, stanfit = fit
    )
  } 
  out
}
