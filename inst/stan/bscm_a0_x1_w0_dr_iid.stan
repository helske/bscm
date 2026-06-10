// Bayesian SCM
// No intercept
// Covariates with constant effects
// Dirichlet prior on omega

data {
#include data/base.stan
#include data/donors.stan
#include data/beta.stan
}
transformed data {
#include transformed_data/omega_dirichlet.stan
}
parameters {
#include parameters/base.stan
#include parameters/omega_dirichlet.stan
#include parameters/beta.stan
}
model {
#include model/base.stan
#include model/omega_dirichlet.stan
#include model/beta.stan
#include model/likelihood_a0_x1_w0_iid.stan
}
generated quantities {
#include generated_quantities/base.stan
#include generated_quantities/y_mean_a0_x1_w0.stan
#include generated_quantities/y_rep_iid.stan
}
