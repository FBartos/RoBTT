functions{
  
#include /chunks/common_functions.stan

}
data {
  // data
  // are individual observation or summary statistics used
  int is_ss;
  // sample sizes
  int<lower=0> N1;
  int<lower=0> N2;
  // individual observations
  vector[is_ss == 0 ? N1 : 0] x1;
  vector[is_ss == 0 ? N2 : 0] x2;
  // summary statistics
  vector[is_ss == 1 ? 2 : 0] mean_i;
  vector[is_ss == 1 ? 2 : 0] sd_i;

  // model type (1 if the parameter is included, 0 if not)
  int is_d;
  int is_r;
  int is_nu;

  // range of the parameters
  vector[is_d  == 1 ? 2 : 0] bounds_d;
  vector[is_r  == 1 ? 2 : 0] bounds_r;
  vector[is_nu == 1 ? 2 : 0] bounds_nu;
  int bounds_type_d[is_d   == 1 ? 2 : 0];
  int bounds_type_r[is_r   == 1 ? 2 : 0];
  int bounds_type_nu[is_nu == 1 ? 2 : 0];

  // prior distribution specification of the parameteres
  real fixed_d[is_d   == 0 ? 1 : 0];
  real fixed_r[is_r   == 0 ? 1 : 0];
  real fixed_nu[is_nu == 0 ? 1 : 0];
  vector[is_d  == 1 ? 3 : 0] prior_parameters_d;
  vector[is_r  == 1 ? 3 : 0] prior_parameters_r;
  vector[is_nu == 1 ? 3 : 0] prior_parameters_nu;
  int prior_type_d;
  int prior_type_r;
  int prior_type_nu;
}
parameters{
  real mu;
  real<lower = 0> scale;
  real<lower = coefs_lb(bounds_type_d[1],  bounds_d[1]),  upper = coefs_ub(bounds_type_d[2],  bounds_d[2])>  d[is_d];
  real<lower = coefs_lb(bounds_type_r[1],  bounds_r[1]),  upper = coefs_ub(bounds_type_r[2],  bounds_r[2])>  r[is_r];
  real<lower = coefs_lb(bounds_type_nu[1], bounds_nu[1]), upper = coefs_ub(bounds_type_nu[2], bounds_nu[2])> nu[is_nu];
}
transformed parameters {
  real pooled_sigma;
  real sigma_i[2];
  real scale_i[2];
  real mu_i[2];
  real nu_t;
  
  
  // compute means and sigmas for each group
  if(is_nu == 1){
    nu_t = nu[1] + 2;
  }else{
    nu_t = fixed_nu[1];
  }
  if(is_r == 1){
    scale_i[1]   = 2 * scale * r[1];
    scale_i[2]   = 2 * scale * (1 - r[1]);
  }else{
    scale_i[1]   = 2 * scale * fixed_r[1];
    scale_i[2]   = 2 * scale * (1 - fixed_r[1]);
  }
  
  sigma_i[1] = scale_i[1] * sqrt(nu_t / (nu_t - 2.0));
  sigma_i[2] = scale_i[2] * sqrt(nu_t / (nu_t - 2.0));
  pooled_sigma = pool_sigma(sigma_i[1], sigma_i[2], N1, N2);
  
  if(is_d == 1){
    mu_i[1] = mu - 0.5 * d[1] * pooled_sigma;
    mu_i[2] = mu + 0.5 * d[1] * pooled_sigma;
  }else{
    mu_i[1] = mu - 0.5 * fixed_d[1] * pooled_sigma;
    mu_i[2] = mu + 0.5 * fixed_d[1] * pooled_sigma;
  }
}
model {
  // default Jeffrey's priors for mu and scale
  target += Jeffreys_mu_lpdf(mu);
  target += Jeffreys_sigma_lpdf(scale);

  // priors on d and r
  if(is_d == 1){
    target += set_prior(d[1], prior_type_d, prior_parameters_d, bounds_type_d, bounds_d);
  }
  if(is_r == 1){
    target += set_prior(r[1], prior_type_r, prior_parameters_r, bounds_type_r, bounds_r);
  }
  if(is_nu == 1){
    target += set_prior(nu[1], prior_type_nu, prior_parameters_nu, bounds_type_nu, bounds_nu);
  }

  // likelihood of the data
  if(is_ss == 0){
    target += student_t_lpdf(x1 | nu_t, mu_i[1], scale_i[1]);
    target += student_t_lpdf(x2 | nu_t, mu_i[2], scale_i[2]);
  }else{
    reject("Fitting models with t likelihood and summary statistics is not possible :(.");
  }
}
