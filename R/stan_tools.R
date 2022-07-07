# prepares stan data object
.stan_data            <- function(data){
  
  if(attr(data, "summary")){
    stan_data <- .stan_data.summary(data[["mean1"]], data[["mean2"]], data[["sd1"]], data[["sd2"]], data[["N1"]], data[["N2"]])
  }else{
    stan_data <- .stan_data.individual(data[["x1"]], data[["x2"]])
  }
  
  return(stan_data)
}
.stan_data.individual <- function(x1, x2){
  return(list(
    x1     = as.array(x1),
    x2     = as.array(x2),
    N1     = length(x1),
    N2     = length(x2),
    is_ss  = 0,
    
    mean_i = numeric(),
    sd_i   = numeric()
  ))
}
.stan_data.summary    <- function(mean1, mean2, sd1, sd2, N1, N2){
  return(list(
    mean_i = c(mean1, mean2),
    sd_i   = c(sd1,   sd2),
    N1     = N1,
    N2     = N2,
    is_ss  = 1,
    
    x1     = numeric(),
    x2     = numeric()
  ))
}

# transforms BayesTools priors into pre-compiled stan code
.stan_distribution            <- function(parameter, prior){
  
  out <- list()
  
  out[[paste0("prior_type_", parameter)]] <- switch(
    prior[["distribution"]],
    "point"      = 0,
    "normal"     = 1,
    "lognormal"  = 2,
    "t"          = 4,
    "gamma"      = 5,
    "invgamma"   = 6,
    "uniform"    = 7,
    "beta"       = 8,
    "exp"        = 9,
  )
  
  if(is.prior.point(prior)){
    
    out[[paste0("is_", parameter)]]          <- 0
    out[[paste0("fixed_", parameter)]]       <- as.array(prior$parameters[["location"]])
    out[[paste0("bounds_", parameter)]]      <- numeric()
    out[[paste0("bounds_type_", parameter)]] <- numeric()
    
    out[[paste0("prior_parameters_", parameter)]] <- numeric()
    
  }else if(is.prior.simple(prior)){
    
    out[[paste0("is_", parameter)]]     <- 1
    out[[paste0("fixed_", parameter)]]  <- numeric()
    
    out[[paste0("bounds_", parameter)]] <- c(
      if(is.infinite(prior$truncation[["lower"]])) 999 else prior$truncation[["lower"]],
      if(is.infinite(prior$truncation[["upper"]])) 999 else prior$truncation[["upper"]]
    )
    
    out[[paste0("bounds_type_", parameter)]] <- c(
      if(is.infinite(prior$truncation[["lower"]])) 0 else 1,
      if(is.infinite(prior$truncation[["upper"]])) 0 else 1
    )
    
    out[[paste0("prior_parameters_", parameter)]] <- .stan_distribution_parameters(prior)
  }else{
    stop("Other prior distributions are not implemented for stan.")
  }
  
  return(out)
}
.stan_distribution_parameters <- function(prior){
  
  # a vector of length three always needs to be passed - filling the redundant values with 999
  prior_parameters <- rep(999, 3)
  
  if(prior[["distribution"]] == "normal"){
    prior_parameters[1] <- prior$parameters[["distribution"]]
    prior_parameters[2] <- prior$parameters[["sd"]]
  }else if(prior[["distribution"]] == "lognormal"){
    prior_parameters[1] <- prior$parameters[["meanlog"]]
    prior_parameters[2] <- prior$parameters[["sdlog"]]
  }else if(prior[["distribution"]] == "t"){
    prior_parameters[1] <- prior$parameters[["df"]]
    prior_parameters[2] <- prior$parameters[["location"]]
    prior_parameters[3] <- prior$parameters[["scale"]]
  }else if(prior[["distribution"]] == "gamma"){
    prior_parameters[1] <- prior$parameters[["shape"]]
    prior_parameters[2] <- prior$parameters[["rate"]]
  }else if(prior[["distribution"]] == "invgamma"){
    prior_parameters[1] <- prior$parameters[["shape"]]
    prior_parameters[2] <- prior$parameters[["scale"]]
  }else if(prior[["distribution"]] == "uniform"){
    prior_parameters[1] <- prior$parameters[["a"]]
    prior_parameters[2] <- prior$parameters[["b"]]
  }else if(prior[["distribution"]] == "beta"){
    prior_parameters[1] <- prior$parameters[["alpha"]]
    prior_parameters[2] <- prior$parameters[["beta"]]
  }else if(prior[["distribution"]] == "exp"){
    prior_parameters[1] <- prior$parameters[["rate"]]
  }
  
  return(prior_parameters)
}