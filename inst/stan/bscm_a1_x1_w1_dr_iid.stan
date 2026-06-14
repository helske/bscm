// Bayesian SCM
// Intercept
// Covariates with varying effects
// Dirichlet prior on omega

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
#include model/likelihood_a1_x1_w1_iid.stan
}
generated quantities {
#include generated_quantities/base.stan
#include generated_quantities/alpha_z.stan
#include generated_quantities/alpha_x.stan
#include generated_quantities/y_mean_a1_x1_w1.stan
#include generated_quantities/y_rep_iid.stan
}
