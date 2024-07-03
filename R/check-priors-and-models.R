.set_priors             <- function(prior_delta, prior_rho, prior_nu, prior_delta_null, prior_rho_null, prior_nu_null, prior_mu, prior_sigma, data = NULL, is_trunc){
  
  priors        <- list()
  priors$delta  <- .set_parameter_priors(prior_delta_null,  prior_delta,  "delta")
  priors$rho    <- .set_parameter_priors(prior_rho_null,    prior_rho,    "rho")
  priors$nu     <- .set_parameter_priors(prior_nu_null,     prior_nu,     "nu")
  
  priors$mu    <- .set_common_prior(prior_mu,    "mu",    data, is_trunc)
  priors$sigma <- .set_common_prior(prior_sigma, "sigma", data, is_trunc)
  
  return(priors)
}
.set_parameter_priors   <- function(priors_null, priors_alt, parameter){
  
  # check that at least one prior is specified (either null or alternative)
  if(parameter != "nu" && (is.null(priors_null) & is.null(priors_alt)))
    stop(paste0("At least one prior needs to be specified for the ", parameter," parameter (either null or alternative)."))
  
  # create an empty list if user didn't specified priors
  if(is.null(priors_null)){
    priors_null <- NULL
  }else{
    # check that the prior is a prior object
    if(is.prior(priors_null)){
      priors_null$is_null <- TRUE
    }else{
      stop(paste0("The null prior distribution for the ", parameter, " was not specified correctly."))
    }
  }
  if(is.null(priors_alt)){
    priors_alt <- NULL
  }else{
    # check that the prior is a prior object
    if(is.prior(priors_alt)){
      priors_alt$is_null <- FALSE
    }else{
      stop(paste0("The alternative prior distribution for the ", parameter, " was not specified correctly."))
    }
  }
  
  # join null and alternative priors
  priors <- list()
  if(!is.null(priors_null)){
    priors$null <- priors_null
  }
  if(!is.null(priors_alt)){
    priors$alt <- priors_alt
  }
  
  
  ### check that the specified prior distributions are valid
  
  if(parameter %in% c("delta")){
    
    # check that the passed priors are supported for the parameter
    if(length(priors) > 0){
      for(p in names(priors)){
        if(!priors[[p]]$distribution %in% c("none", "normal", "lognormal", "t", "gamma", "invgamma", "point", "uniform", "beta", "exp"))
          stop(paste0(priors[[p]]$distribution," prior distribution is not supported for the ", parameter," parameter. See '?prior' for further information."))
        if(priors[[p]]$distribution == "none"){
          temp_is_null        <- priors[[p]]$is_null
          priors[[p]]         <- prior(distribution = "spike", parameters = list(location = 0), prior_weights = priors[[p]][["prior_weights"]])
          priors[[p]]$is_null <- temp_is_null
        }
      }
    }
    
    
  }else if(parameter == "rho"){
    
    # check that the passed priors are supported for the r parameter
    if(length(priors) > 0){
      for(p in names(priors)){
        if(!priors[[p]]$distribution %in% c("none", "normal", "lognormal", "t", "gamma", "invgamma", "point", "uniform", "beta", "exp"))
          stop(paste0(priors[[p]]$distribution," prior distribution is not supported for the tau parameter. See '?prior' for further information."))
        if(priors[[p]]$distribution == "none"){
          temp_is_null        <- priors[[p]]$is_null
          priors[[p]]         <- prior(distribution = "spike", parameters = list(location = 0.5), prior_weights = priors[[p]][["prior_weights"]])
          priors[[p]]$is_null <- temp_is_null
        }else if(priors[[p]]$distribution == "point"){
          if(priors[[p]]$parameters$location <= 0 | priors[[p]]$parameters$location >= 1){
            stop(paste0("The location of a point prior distribution for ", parameter, " parameter must be larger than 0 and lower than 1. See '?prior' for further information."))
          }
        }else if(priors[[p]]$distribution == "uniform"){
          if(priors[[p]]$parameters$b < 0 ){
            stop(paste0("The uniform prior distribution for ", parameter, " parameter cannot be defined on negative numbers. See '?prior' for further information."))
          }
          if(priors[[p]]$parameters$a > 1){
            stop(paste0("The uniform prior distribution for ", parameter, " parameter cannot be defined on numbers larger than 1."))
          }
        }else{
          if(priors[[p]]$truncation$lower < 0){
            priors[[p]]$truncation$lower <- 0
            warning(paste0("The range of a prior distribution for ", parameter, " parameter cannot contain values lower than 0. The lower truncation point was set to 0. See '?prior' for further information."))
          }
          if(priors[[p]]$truncation$upper > 1){
            priors[[p]]$truncation$upper <- 1
            warning(paste0("The range of a prior distribution for ", parameter, " parameter cannot contain values larger than 1. The upper truncation point was set to 1. See '?prior' for further information."))
          }
        }
      }
    }
    
    
  }else if(parameter == "nu"){
    
    # check that the passed priors are supported for the nu parameter
    if(length(priors) > 0){
      for(p in names(priors)){
        if(!priors[[p]]$distribution %in% c("none", "normal", "lognormal", "t", "gamma", "invgamma", "point", "uniform", "beta", "exp"))
          stop(paste0(priors[[p]]$distribution," prior distribution is not supported for the tau parameter. See '?prior' for further information."))
        if(priors[[p]]$distribution == "none")
          next
        if(priors[[p]]$distribution == "point"){
          if(is.infinite(priors[[p]]$parameters$location)){
            temp_is_null        <- priors[[p]]$is_null
            priors[[p]]         <- prior_none(prior_weights = priors[[p]][["prior_weights"]])
            priors[[p]]$is_null <- temp_is_null
          }else if(priors[[p]]$parameters$location <= 0){
            stop(paste0("The location of a point prior distribution for ", parameter, " parameter must be larger than 0 and lower than 1. See '?prior' for further information."))
          }
        }else if(priors[[p]]$distribution == "uniform"){
          if(priors[[p]]$parameters$b < 0){
            stop(paste0("The uniform prior distribution for ", parameter, " parameter cannot be defined on negative numbers. See '?prior' for further information."))
          }
        }else{
          if(priors[[p]]$truncation$lower < 0){
            priors[[p]]$truncation$lower <- 0
            warning(paste0("The range of a prior distribution for ", parameter, " parameter cannot contain values lower than 0. The lower truncation point was set to 0. See '?prior' for further information."))
          }
        }
      }
    }
  }
  
  
  return(priors)
}
.set_common_prior       <- function(prior, parameter, data = NULL, is_trunc){
  
  # set default priors if not specified
  if(is.null(prior)){
    if(is_trunc){
      if(is.null(data)) {
        prior <- switch(
          parameter,
          "mu"    = prior("cauchy", parameters = list(location = 0, scale = 1)),
          "sigma" = prior("rate", parameters = list(rate = 1))
        )
      } else {
        
        if (attr(data, "summary")) {
          mean_x <- (data[["mean1"]] + data[["mean2"]]) / 2
          var_x  <- (data[["sd1"]]^2 + data[["sd2"]]^2) / 2
        } else {
          mean_x <- (mean(data[["x1"]]) + mean(data[["x2"]])) / 2
          var_x  <- (var(data[["x1"]])  + var(data[["x2"]]))  / 2
        }
        
        prior <- switch(
          parameter,
          "mu"    = prior("cauchy", parameters = list(location = mean_x, scale = sqrt(var_x))),
          "sigma" = prior("exp",    parameters = list(rate = 1/sqrt(var_x)))
        )
      }
    }else{
      prior <- switch(
        parameter,
        "mu"    = "Jeffreys_mu",
        "sigma" = "Jeffreys_sigma"
      )
    }
    return(prior)
  }
  
  ### check that the specified prior distributions are valid
  if(parameter == "mu"){
    
    if(!prior$distribution %in% c("normal", "lognormal", "t", "gamma", "invgamma", "uniform", "beta", "exp"))
      stop(paste0(prior$distribution," prior distribution is not supported for the ", parameter," parameter. See '?prior' for further information."))
    
  }else if(parameter %in% "sigma"){
    
    if(!prior$distribution %in% c("normal", "lognormal", "t", "gamma", "invgamma", "uniform", "beta", "exp"))
      stop(paste0(prior$distribution," prior distribution is not supported for the ", parameter, " parameter. See '?prior' for further information."))
    if(prior$distribution == "uniform"){
      if(prior$parameters$a < 0 ){
        stop(paste0("The uniform prior distribution for ", parameter, " parameter cannot be defined on negative numbers. See '?prior' for further information."))
      }
    }else{
      if(prior$truncation$lower < 0){
        prior$truncation$lower <- 0
        warning(paste0("The range of a prior distribution for ", parameter, " parameter cannot contain values lower than 0. The lower truncation point was set to 0. See '?prior' for further information."))
      }
    }
  }
  
  return(prior)
}
.get_models             <- function(priors){
  
  # create models according to the set priors
  models <- NULL
  for(delta in priors[["delta"]]){
    for(rho in priors[["rho"]]){
      for(nu in priors[["nu"]]){
        models <- c(models, list(.create_model(delta, rho, nu, NULL, NULL, delta[["prior_weights"]] * rho[["prior_weights"]] * nu[["prior_weights"]], priors[["mu"]], priors[["sigma"]])))
      }
    }
  }
  
  return(models)
}
.create_model           <- function(prior_delta, prior_rho, prior_nu, prior_p, likelihood, prior_weights, prior_mu, prior_sigma){
  
  priors <- list()
  
  priors$delta <- prior_delta
  priors$rho   <- prior_rho
  priors$nu    <- prior_nu

  priors$mu    <- prior_mu
  priors$sigma <- prior_sigma
  
  # possibly simplify t to normal
  if(prior_nu$distribution == "none"){
    likelihood <- "normal"
  }else{
    likelihood <- "t"
  }
  
  model <- list(
    priors            = priors,
    prior_weights     = prior_weights,
    prior_weights_set = prior_weights,
    likelihood        = likelihood
  )
  class(model) <- "RoBTT.model"
  
  return(model)
}
