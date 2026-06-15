// Bayesian SCM
// No intercept
// Covariates with constant effects

functions {
#include functions/omega_logistic_normal.stan
}
data {
#include data/base.stan
#include data/beta.stan
}
transformed data {
#include transformed_data/base.stan
#include transformed_data/omega.stan
}
parameters {
#include parameters/base.stan
#include parameters/omega.stan
#include parameters/beta.stan
}
transformed parameters {
#include transformed_parameters/omega.stan
}
model {
#include model/base.stan
#include model/omega.stan
#include model/beta.stan
#include model/likelihood_a0_x1_w0_iid.stan
}
generated quantities {
#include generated_quantities/base.stan
#include generated_quantities/y_mean_a0_x1_w0.stan
#include generated_quantities/y_rep_iid.stan
}
