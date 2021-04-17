functions{
  // default Jeffrey priors for sigma and mu
  real Jeffreys_mu_lpdf(real mu){
    return 0;
  }
  real Jeffreys_sigma_lpdf(real sigma){
    return log(1/sigma);
  }

  // function for computation of pooled standard deviation
  real pool_sigma(real sigma1, real sigma2, real N1, real N2){
    return sqrt( (pow(sigma1, 2) * (N1 - 1) + pow(sigma2, 2) * (N2 - 1)) / (N1 + N2 - 2) );
  }

  // function for setting parameter bounds
  real coefs_lb(int type, real bound) {
    real lb;
    if (type == 0)
      lb = negative_infinity();
    else
      lb = bound;
    return lb;
  }
  real coefs_ub(int type, real bound) {
    real lb;
    if (type == 0)
      lb = positive_infinity();
    else
      lb = bound;
    return lb;
  }

  // function for setting prior distributions on the parameters
  // type 1 = normal
  // type 2 = lognormal
  // type 3 = cauchy
  // type 4 = t
  // type 5 = gamma
  // type 6 = inverse-gamma
  // type 7 = uniform
  // type 8 = beta
  real set_prior(real parameter, int prior_type, vector prior_parameters, int[] bounds_type, vector bounds){
    real ll;

    if(prior_type == 1){
      ll = normal_lpdf(parameter | prior_parameters[1], prior_parameters[2]);
      if(bounds_type[1] != 0 && bounds_type[2] != 0){
        ll -= log_diff_exp(
          normal_lcdf(bounds[2] | prior_parameters[1], prior_parameters[2]),
          normal_lcdf(bounds[1] | prior_parameters[1], prior_parameters[2])
        );
      }else if(bounds_type[1] != 0){
        ll -= normal_lccdf(bounds[1] | prior_parameters[1], prior_parameters[2]);
      }else if(bounds_type[2] != 0){
        ll -= normal_lcdf(bounds[2]  | prior_parameters[1], prior_parameters[2]);
      }
    }else if(prior_type == 2){
      ll = lognormal_lpdf(parameter | prior_parameters[1], prior_parameters[2]);
      if(bounds_type[1] != 0 && bounds_type[2] != 0){
        ll -= log_diff_exp(
          lognormal_lcdf(bounds[2] | prior_parameters[1], prior_parameters[2]),
          lognormal_lcdf(bounds[1] | prior_parameters[1], prior_parameters[2])
        );
      }else if(bounds_type[1] != 0){
        ll -= lognormal_lccdf(bounds[1] | prior_parameters[1], prior_parameters[2]);
      }else if(bounds_type[2] != 0){
        ll -= lognormal_lcdf(bounds[2]  | prior_parameters[1], prior_parameters[2]);
      }
    }else if(prior_type == 3){
      ll = cauchy_lpdf(parameter | prior_parameters[1], prior_parameters[2]);
      if(bounds_type[1] != 0 && bounds_type[2] != 0){
        ll -= log_diff_exp(
          cauchy_lcdf(bounds[2] | prior_parameters[1], prior_parameters[2]),
          cauchy_lcdf(bounds[1] | prior_parameters[1], prior_parameters[2])
        );
      }else if(bounds_type[1] != 0){
        ll -= cauchy_lccdf(bounds[1] | prior_parameters[1], prior_parameters[2]);
      }else if(bounds_type[2] != 0){
        ll -= cauchy_lcdf(bounds[2]  | prior_parameters[1], prior_parameters[2]);
      }
    }else if(prior_type == 4){
      ll = student_t_lpdf(parameter | prior_parameters[1], prior_parameters[2], prior_parameters[3]);
      if(bounds_type[1] != 0 && bounds_type[2] != 0){
        ll -= log_diff_exp(
          student_t_lcdf(bounds[2] | prior_parameters[1], prior_parameters[2], prior_parameters[3]),
          student_t_lcdf(bounds[1] | prior_parameters[1], prior_parameters[2], prior_parameters[3])
        );
      }else if(bounds_type[1] != 0){
        ll -= student_t_lccdf(bounds[1] | prior_parameters[1], prior_parameters[2], prior_parameters[3]);
      }else if(bounds_type[2] != 0){
        ll -= student_t_lcdf(bounds[2]  | prior_parameters[1], prior_parameters[2], prior_parameters[3]);
      }
    }else if(prior_type == 5){
      ll = gamma_lpdf(parameter | prior_parameters[1], prior_parameters[2]);
      if(bounds_type[1] != 0 && bounds_type[2] != 0){
        ll -= log_diff_exp(
          gamma_lcdf(bounds[2] | prior_parameters[1], prior_parameters[2]),
          gamma_lcdf(bounds[1] | prior_parameters[1], prior_parameters[2])
        );
      }else if(bounds_type[1] != 0){
        ll -= gamma_lccdf(bounds[1] | prior_parameters[1], prior_parameters[2]);
      }else if(bounds_type[2] != 0){
        ll -= gamma_lcdf(bounds[2]  | prior_parameters[1], prior_parameters[2]);
      }
    }else if(prior_type == 6){
      ll = inv_gamma_lpdf(parameter | prior_parameters[1], prior_parameters[2]);
      if(bounds_type[1] != 0 && bounds_type[2] != 0){
        ll -= log_diff_exp(
          inv_gamma_lcdf(bounds[2] | prior_parameters[1], prior_parameters[2]),
          inv_gamma_lcdf(bounds[1] | prior_parameters[1], prior_parameters[2])
        );
      }else if(bounds_type[1] != 0){
        ll -= inv_gamma_lccdf(bounds[1] | prior_parameters[1], prior_parameters[2]);
      }else if(bounds_type[2] != 0){
        ll -= inv_gamma_lcdf(bounds[2]  | prior_parameters[1], prior_parameters[2]);
      }
    }else if(prior_type == 7){
      ll = uniform_lpdf(parameter | prior_parameters[1], prior_parameters[2]);
    }else if(prior_type == 8){
      ll = beta_lpdf(parameter | prior_parameters[1], prior_parameters[2]);
      if(bounds_type[1] != 0 && bounds_type[2] != 0){
        ll -= log_diff_exp(
          beta_lcdf(bounds[2] | prior_parameters[1], prior_parameters[2]),
          beta_lcdf(bounds[1] | prior_parameters[1], prior_parameters[2])
        );
      }else if(bounds_type[1] != 0){
        ll -= beta_lccdf(bounds[1] | prior_parameters[1], prior_parameters[2]);
      }else if(bounds_type[2] != 0){
        ll -= beta_lcdf(bounds[2]  | prior_parameters[1], prior_parameters[2]);
      }
    }

    return ll;
  }

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

  // range of the parameters
  vector[2] bounds_d;
  vector[2] bounds_r;
  int bounds_type_d[2];
  int bounds_type_r[2];

  // prior distribution specification of the parameteres
  vector[3] prior_parameters_d;
  vector[3] prior_parameters_r;
  int prior_type_d;
  int prior_type_r;
}
parameters{
  real mu;
  real<lower = 0> sigma;
  real<lower = coefs_lb(bounds_type_d[1], bounds_d[1]), upper = coefs_ub(bounds_type_d[2], bounds_d[2])> d[is_d];
  real<lower = coefs_lb(bounds_type_r[1], bounds_r[1]), upper = coefs_ub(bounds_type_r[2], bounds_r[2])> r[is_r];
}
transformed parameters {
  real pooled_sigma;
  real sigma_i[2];
  real mu_i[2];

  // compute means and sigmas for each group
  if(is_r == 1){
    sigma_i[1]   = sigma * (2 * r[1]);
    sigma_i[2]   = sigma * (2 * (1 - r[1]));
    pooled_sigma = pool_sigma(sigma_i[1], sigma_i[2], N1, N2);
  }else{
    sigma_i[1]   = sigma;
    sigma_i[2]   = sigma;
    pooled_sigma = sigma;
  }
  if(is_d == 1){
    mu_i[1] = mu - 0.5 * d[1] * pooled_sigma;
    mu_i[2] = mu + 0.5 * d[1] * pooled_sigma;
  }else{
    mu_i[1] = mu;
    mu_i[2] = mu;
  }
}
model {
  // default Jeffrey's priors for mu and sigma
  target += Jeffreys_mu_lpdf(mu);
  target += Jeffreys_sigma_lpdf(sigma);

  // priors on d and r
  if(is_d == 1){
    target += set_prior(d[1], prior_type_d, prior_parameters_d, bounds_type_d, bounds_d);
  }
  if(is_r == 1){
    target += set_prior(r[1], prior_type_r, prior_parameters_r, bounds_type_r, bounds_r);
  }

  // likelihood of the data
  if(is_ss == 0){
    target += normal_lpdf(x1 | mu_i[1], sigma_i[1]);
    target += normal_lpdf(x2 | mu_i[2], sigma_i[2]);
  }else{
    target += normal_lpdf(mean_i[1] | mu_i[1], sigma_i[1] / sqrt(N1));
    target += normal_lpdf(mean_i[2] | mu_i[2], sigma_i[2] / sqrt(N2));
    target += gamma_lpdf((N1 - 1) * pow(sd_i[1], 2) | 0.5 * (N1 - 1), 1 / (2 * pow(sigma_i[1], 2)));
    target += gamma_lpdf((N2 - 1) * pow(sd_i[2], 2) | 0.5 * (N2 - 1), 1 / (2 * pow(sigma_i[2], 2)));
  }
}
