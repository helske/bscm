// Bayesian SCM
// No intercept
// Covariates with constant effects
// Dirichlet prior on omega
// AR(1) residuals

data {
#include data/base.stan
#include data/donors.stan
#include data/beta.stan
}
transformed data {
#include transformed_data/base.stan
#include transformed_data/omega_dirichlet.stan
}
parameters {
#include parameters/base.stan
#include parameters/omega_dirichlet.stan
#include parameters/beta.stan
#include parameters/rho.stan
}
model {
#include model/base.stan
#include model/omega_dirichlet.stan
#include model/beta.stan
#include model/rho.stan
#include model/likelihood_a0_x1_w0_ar1.stan
}
generated quantities {
#include generated_quantities/base.stan
#include generated_quantities/gq_ar1_a0_x1_w0.stan
}
