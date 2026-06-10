// Intercept-only model, no donors
data {
#include data/base.stan
#include data/alpha.stan
}
parameters {
#include parameters/base.stan
#include parameters/a.stan
}
transformed parameters {
  vector[N] alpha = a;
}
model {
#include model/base.stan
#include model/a.stan
#include model/likelihood_a1_x0_w0_no_iid.stan
}
generated quantities {
#include generated_quantities/base.stan
#include generated_quantities/y_mean_a1_x0_w0_no.stan
#include generated_quantities/y_rep_iid.stan
}
