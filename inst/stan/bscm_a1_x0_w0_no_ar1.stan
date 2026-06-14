// Intercept-only model, no donors
// AR(1) residuals
data {
#include data/base.stan
#include data/alpha.stan
}
parameters {
#include parameters/base.stan
#include parameters/a.stan
#include parameters/rho.stan
}
model {
#include model/base.stan
#include model/a.stan
#include model/rho.stan
#include model/likelihood_a1_x0_w0_no_ar1.stan
}
generated quantities {
#include generated_quantities/base.stan
  vector[N] alpha = a;
#include generated_quantities/gq_ar1_a1_x0_w0_no.stan
}
