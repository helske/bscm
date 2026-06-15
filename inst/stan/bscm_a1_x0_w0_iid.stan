// Bayesian SCM
// Intercept
// No covariates

functions {
#include functions/omega_logistic_normal.stan
}
data {
#include data/base.stan
#include data/alpha.stan
}
transformed data {
#include transformed_data/base.stan
#include transformed_data/omega.stan
#include transformed_data/z_means.stan
}
parameters {
#include parameters/base.stan
#include parameters/omega.stan
#include parameters/a.stan
}
transformed parameters {
#include transformed_parameters/omega.stan
}
model {
#include model/base.stan
#include model/omega.stan
#include model/likelihood_a1_x0_w0_iid.stan
}
generated quantities {
#include generated_quantities/base.stan
#include generated_quantities/alpha_z.stan
#include generated_quantities/y_mean_a1_x0_w0.stan
#include generated_quantities/y_rep_iid.stan
}
