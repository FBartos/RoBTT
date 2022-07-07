.set_priors             <- function(prior_delta, prior_rho, prior_nu, prior_delta_null, prior_rho_null, prior_nu_null){
  
  priors        <- list()
  priors$delta  <- .set_parameter_priors(prior_delta_null,  prior_delta,  "delta")
  priors$rho    <- .set_parameter_priors(prior_rho_null,    prior_rho,    "rho")
  priors$nu     <- .set_parameter_priors(prior_nu_null,     prior_nu,     "nu")
  
  return(priors)
}
.set_parameter_priors   <- function(priors_null, priors_alt, parameter){
  
  # check that at least one prior is specified (either null or alternative)
  if(is.null(priors_null) & is.null(priors_alt))
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
      for(i in 1:length(priors)){
        if(!priors[[i]]$distribution %in% c("normal", "lognormal", "t", "gamma", "invgamma", "point", "uniform", "beta", "exp"))
          stop(paste0(priors[[i]]$distribution," prior distribution is not supported for the ", parameter," parameter. See '?prior' for further information."))
      }
    }
    
    
  }else if(parameter == "rho"){
    
    # check that the passed priors are supported for the r parameter
    if(length(priors) > 0){
      for(i in 1:length(priors)){
        if(!priors[[i]]$distribution %in% c("normal", "lognormal", "t", "gamma", "invgamma", "point", "uniform", "beta", "exp"))stop(paste0(priors[[i]]$distribution," prior distribution is not supported for the tau parameter. See '?prior' for further information."))
        if(priors[[i]]$distribution == "point"){
          if(priors[[i]]$parameters$location <= 0 | priors[[i]]$parameters$location >= 1){
            stop(paste0("The location of a point prior distribution for ", parameter, " parameter must be larger than 0 and lower than 1. See '?prior' for further information."))
          }
        }else if(priors[[i]]$distribution == "uniform"){
          if(priors[[i]]$parameters$b < 0 ){
            stop(paste0("The uniform prior distribution for ", parameter, " parameter cannot be defined on negative numbers. See '?prior' for further information."))
          }
          if(priors[[i]]$parameters$a > 1){
            stop(paste0("The uniform prior distribution for ", parameter, " parameter cannot be defined on numbers larger than 1."))
          }
        }else{
          if(priors[[i]]$truncation$lower < 0){
            priors[[i]]$truncation$lower <- 0
            warning(paste0("The range of a prior distribution for ", parameter, " parameter cannot contain values lower than 0. The lower truncation point was set to 0. See '?prior' for further information."))
          }
          if(priors[[i]]$truncation$upper > 1){
            priors[[i]]$truncation$upper <- 1
            warning(paste0("The range of a prior distribution for ", parameter, " parameter cannot contain values larger than 1. The upper truncation point was set to 1. See '?prior' for further information."))
          }
        }
      }
    }
    
    
  }else if(parameter == "nu"){
    
    # check that the passed priors are supported for the nu parameter
    if(length(priors) > 0){
      for(i in 1:length(priors)){
        if(!priors[[i]]$distribution %in% c("normal", "lognormal", "t", "gamma", "invgamma", "point", "uniform", "beta", "exp"))stop(paste0(priors[[i]]$distribution," prior distribution is not supported for the tau parameter. See '?prior' for further information."))
        if(priors[[i]]$distribution == "point"){
          if(priors[[i]]$parameters$location <= 0){
            stop(paste0("The location of a point prior distribution for ", parameter, " parameter must be larger than 0 and lower than 1. See '?prior' for further information."))
          }
        }else if(priors[[i]]$distribution == "uniform"){
          if(priors[[i]]$parameters$b < 0 ){
            stop(paste0("The uniform prior distribution for ", parameter, " parameter cannot be defined on negative numbers. See '?prior' for further information."))
          }
        }else{
          if(priors[[i]]$truncation$lower < 0){
            priors[[i]]$truncation$lower <- 0
            warning(paste0("The range of a prior distribution for ", parameter, " parameter cannot contain values lower than 0. The lower truncation point was set to 0. See '?prior' for further information."))
          }
        }
      }
    }
  }
  
  
  return(priors)
}
.get_models             <- function(priors, likelihoods){
  
  # create models according to the set priors
  models <- NULL
  for(delta in priors[["delta"]]){
    for(rho in priors[["rho"]]){
      for(likelihood in likelihoods){
        if(likelihood == "t"){
          for(nu in priors$nu){
            models <- c(models, list(.create_model(delta, rho, nu, NULL, likelihood, delta[["prior_weights"]] * rho[["prior_weights"]] * nu[["prior_weights"]])))
          }
        }else{
          models <- c(models, list(.create_model(delta, rho, NULL, NULL, likelihood, delta[["prior_weights"]] * rho[["prior_weights"]])))
        } 
      }
    }
  }
  
  
  return(models)
}
.create_model           <- function(prior_delta, prior_rho, prior_nu, prior_p, likelihood, prior_weights){
  
  priors <- list()
  
  priors$delta <- prior_delta
  priors$rho   <- prior_rho
  priors$nu    <- prior_nu
  
  # possibly simplify t to normal
  if(likelihood == "t" && prior_nu$distribution == "point" && prior_nu$parameters[["location"]] == Inf){
    likelihood <- "normal"
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
