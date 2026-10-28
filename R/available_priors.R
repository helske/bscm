#' Available prior distributions
#'
#' Lists the prior distributions which can be used in the `priors` argument of
#' [bscm()]. See [bscm_prior()] for which parameters of the model each prior
#' can be used with.
#'
#' The following priors follow the function definitions used in Stan:
#'
#' * `beta_pr(shape1, shape2)`: Beta distribution
#' * `cauchy_pr(location, scale)`: Cauchy distribution
#' * `exponential_pr(rate)`: Exponential distribution
#' * `gamma_pr(shape, rate)`: Gamma distribution
#' * `inv_gamma_pr(shape, scale)`: Inverse gamma distribution
#' * `laplace_pr(location, scale)`: Laplace (double exponential) distribution
#' * `normal_pr(location, scale)`: Normal distribution
#' * `student_pr(df, location, scale)`: Student's t distribution
#' 
#' In addition, following priors are available:
#' * `half_normal_pr(scale)`: Half-normal distribution. This corresponds to 
#'   `normal_pr(0, scale)` prior for parameter with lower bound of zero.
#' * `sqrt_inv_gamma_pr(shape, scale)`: Square-root inverse gamma
#'   distribution defined as
#'   \eqn{x^2 \sim \text{Inv-Gamma}(\text{shape}, \text{scale})}, with the
#'   inverse gamma distribution as in `inv_gamma_pr()`.
#' * `scaled_inv_chi_pr(df, scale)`: Scaled inverse chi distribution defined
#'   as 
#'   \eqn{x^2 \sim \text{Inv-Gamma}(\text{df} / 2, \text{scale}^2 / 2)}, i.e.
#'   `sqrt_inv_gamma_pr(df / 2, scale^2 / 2)`.
#' 
#' For the priors on donor weights \eqn{\omega}, following options are 
#' available:
#' * `dirichlet_pr(concentration)`: Symmetric Dirichlet distribution with
#'   concentration \eqn{\zeta}, i.e., 
#'   \eqn{\omega \sim \text{Dirichlet}(\zeta, \ldots, \zeta)}.
#' * `logistic_normal_pr(concentration)`: Logistic normal distribution where 
#'   concentration \eqn{\zeta} is the precision of the the normal distribution, 
#'   so that \eqn{\omega = \text{softmax}(\eta/\sqrt(\zeta))}, where
#'   \eqn{\eta \sim N(0, I)} constrained to sum to zero and scaled to have 
#'   marginal standard deviation of one.
#' * `logistic_laplace_pr(concentration)`: As the logistic normal prior, but
#'   with the zero-sum vector \eqn{\eta} constructed from independent (scaled) 
#'   Laplace-distributed variables.
#' 
#' For all of the three priors above, small values of concentration corresponds 
#' to prior favoring sparse weight vectors, while larger concentration 
#' correspond to more diffuse simplex vectors. Specifically for Dirichlet prior,
#' \eqn{\zeta < 1} concentrate weight on 
#' few donors, while \eqn{\zeta = 1} is uniform prior over the simplex, and
#' \eqn{\zeta > 1} pulls weights toward the center of the simplex.
#' 
#' The default prior on the concentration \eqn{\zeta} is Inverse-Gamma with
#' shape \eqn{\sqrt(J)} and scale \eqn{\sqrt(J) / (2\log(1 + J / 5))}, where 
#' \eqn{J} is the number of donors.
#'
#' @return A `data.frame` with columns `distribution` (name of the
#'   distribution), `constructor` (function defining the prior),
#'   `hyperparameters` (arguments of the constructor), and `support`
#'   (support of the distribution).
#' @seealso [bscm_prior()], [bscm()]
#' @references Stan Development Team. 2026 Stan Function Reference, 2.40. 
#' https://mc-stan.orgStan https://mc-stan.org/docs/functions-reference
#' @export
#' @examples
#' available_priors()
available_priors <- function() {
  data.frame(
    distribution = c(
      "beta",
      "cauchy",
      "dirichlet",
      "exponential",
      "gamma",
      "half_normal",
      "inv_gamma",
      "laplace",
      "logistic_laplace",
      "logistic_normal",
      "normal",
      "scaled_inv_chi",
      "sqrt_inv_gamma",
      "student_t"
    ),
    constructor = c(
      "beta_pr()",
      "cauchy_pr()",
      "dirichlet_pr()",
      "exponential_pr()",
      "gamma_pr()",
      "half_normal_pr()",
      "inv_gamma_pr()",
      "laplace_pr()",
      "logistic_laplace_pr()",
      "logistic_normal_pr()",
      "normal_pr()",
      "scaled_inv_chi_pr()",
      "sqrt_inv_gamma_pr()",
      "student_pr()"
    ),
    hyperparameters = c(
      "shape1, shape2",
      "location, scale",
      "concentration",
      "rate",
      "shape, rate",
      "scale",
      "shape, scale",
      "location, scale",
      "scale",
      "scale",
      "location, scale",
      "df, scale",
      "shape, scale",
      "df, location, scale"
    ),
    support = c(
      "(0, 1)",
      "real",
      "simplex",
      "positive",
      "positive",
      "positive",
      "positive",
      "real",
      "simplex",
      "simplex",
      "real",
      "positive",
      "positive",
      "real"
    )
  )
}
