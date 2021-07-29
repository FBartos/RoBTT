 // default Jeffrey priors for sigma and mu
  real Jeffreys_mu_lpdf(real mu){
    return 0;
  }
  real Jeffreys_sigma_lpdf(real sigma2){
    return log(1/sigma2);
  }

  // function for computation of pooled standard deviation
  real pool_sigma(real sigma1, real sigma2, real N1, real N2){
    return sqrt( (pow(sigma1, 2) * N1 + pow(sigma2, 2) * N2) / (N1 + N2) );
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
    }else if(prior_type == 9){
      ll = exponential_lpdf(parameter | prior_parameters[1]);
      if(bounds_type[1] != 0 && bounds_type[2] != 0){
        ll -= log_diff_exp(
          exponential_lcdf(bounds[2] | prior_parameters[1]),
          exponential_lcdf(bounds[1] | prior_parameters[1])
        );
      }else if(bounds_type[1] != 0){
        ll -= exponential_lccdf(bounds[1] | prior_parameters[1]);
      }else if(bounds_type[2] != 0){
        ll -= exponential_lcdf(bounds[2]  | prior_parameters[1]);
      }
    }

    return ll;
  }