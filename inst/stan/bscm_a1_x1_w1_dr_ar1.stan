// Bayesian SCM
// Intercept
// Covariates with varying effects
// Dirichlet prior on omega
// AR(1) residuals

data {
#include data/base.stan
#include data/donors.stan
#include data/alpha.stan
#include data/beta.stan
#include data/gamma.stan
}
transformed data {
#include transformed_data/base.stan
#include transformed_data/omega_dirichlet.stan
#include transformed_data/z_means.stan
#include transformed_data/x_means.stan
#include transformed_data/gamma.stan
}
parameters {
#include parameters/base.stan
#include parameters/omega_dirichlet.stan
#include parameters/a.stan
#include parameters/beta.stan
#include parameters/gamma.stan
#include parameters/rho.stan
}
transformed parameters {
#include transformed_parameters/gamma.stan
}
model {
#include model/base.stan
#include model/omega_dirichlet.stan
#include model/a.stan
#include model/beta.stan
#include model/gamma.stan
#include model/rho.stan
#include model/likelihood_a1_x1_w1_ar1.stan
}
generated quantities {
#include generated_quantities/base.stan
#include generated_quantities/alpha_z.stan
#include generated_quantities/alpha_x.stan
#include generated_quantities/gq_ar1_a1_x1_w1.stan
}
