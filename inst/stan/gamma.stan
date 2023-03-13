#include /include/common_functions.stan

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

  // range of the parameters
  vector[is_d == 1 ? 2 : 0] bounds_d;
  vector[is_r == 1 ? 2 : 0] bounds_r;
  int bounds_type_d[is_d == 1 ? 2 : 0];
  int bounds_type_r[is_r == 1 ? 2 : 0];

  // prior distribution specification of the parameteres
  vector[is_d == 0 ? 1 : 0] fixed_d;
  vector[is_r == 0 ? 1 : 0] fixed_r;
  vector[is_d == 1 ? 3 : 0] prior_parameters_d;
  vector[is_r == 1 ? 3 : 0] prior_parameters_r;
  int prior_type_d;
  int prior_type_r;
}
parameters{
  real<lower = 0> mu;
  real<lower = 0> sigma2;
  real<lower = coefs_lb(bounds_type_d[1], bounds_d[1]), upper = coefs_ub(bounds_type_d[2], bounds_d[2])> delta[is_d];
  real<lower = coefs_lb(bounds_type_r[1], bounds_r[1]), upper = coefs_ub(bounds_type_r[2], bounds_r[2])> rho[is_r];
}
transformed parameters {
  real pooled_sigma;
  vector[2] mu_i;
  vector[2] sigma_i;
  vector[2] alpha_i;
  vector[2] beta_i;

  // compute means and sigmas for each group
  if(is_r == 1){
    sigma_i[1]   = sqrt( 1 / (2 * 1/sigma2 * rho[1]       ) );
    sigma_i[2]   = sqrt( 1 / (2 * 1/sigma2 * (1 - rho[1]) ) );
    pooled_sigma = pool_sigma(sigma_i[1], sigma_i[2], N1, N2);
  }else{
    sigma_i[1]   = sqrt( 1 / (2 * 1/sigma2 * fixed_r[1]       ) );
    sigma_i[2]   = sqrt( 1 / (2 * 1/sigma2 * (1 - fixed_r[1]) ) );
    pooled_sigma = pool_sigma(sigma_i[1], sigma_i[2], N1, N2);
  }
  if(is_d == 1){
    mu_i[1] = mu - 0.5 * delta[1] * pooled_sigma;
    mu_i[2] = mu + 0.5 * delta[1] * pooled_sigma;
  }else{
    mu_i[1] = mu - 0.5 * fixed_d[1] * pooled_sigma;
    mu_i[2] = mu + 0.5 * fixed_d[1] * pooled_sigma;
  }

  for(i in 1:2){
    alpha_i[i] = pow(mu_i[i], 2) / pow(sigma_i[i], 2);
    beta_i[i]  = mu_i[i] / pow(sigma_i[i], 2);
  }
}
model {
  // default Jeffrey's priors for mu and sigma
  target += Jeffreys_mu_lpdf(mu);
  target += Jeffreys_sigma_lpdf(sigma2);

  // priors on d and r
  if(is_d == 1){
    target += set_prior(delta[1], prior_type_d, prior_parameters_d, bounds_type_d, bounds_d);
  }
  if(is_r == 1){
    target += set_prior(rho[1], prior_type_r, prior_parameters_r, bounds_type_r, bounds_r);
  }

  // likelihood of the data
  if(is_ss == 0){
    target += gamma_lpdf(x1 | alpha_i[1], beta_i[1]);
    target += gamma_lpdf(x2 | alpha_i[2], beta_i[2]);
  }else{
    reject("Fitting models with t likelihood and summary statistics is not possible :(.");
  }
}
