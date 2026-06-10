// Bayesian SCM
// No intercept
// Covariates with constant effects
// Logistic normal prior on omega
// AR(1) residuals
functions {
#include functions/omega_logistic_normal.stan
}
data {
#include data/base.stan
#include data/donors.stan
#include data/beta.stan
}
transformed data {
#include transformed_data/omega_logistic_normal.stan
}
parameters {
#include parameters/base.stan
#include parameters/omega_logistic_normal.stan
#include parameters/beta.stan
#include parameters/rho.stan
}
transformed parameters {
#include transformed_parameters/omega_logistic_normal.stan
}
model {
#include model/base.stan
#include model/omega_logistic_normal.stan
#include model/beta.stan
#include model/rho.stan
#include model/likelihood_a0_x1_w0_ar1.stan
}
generated quantities {
#include generated_quantities/base.stan
#include generated_quantities/gq_ar1_a0_x1_w0.stan
}
