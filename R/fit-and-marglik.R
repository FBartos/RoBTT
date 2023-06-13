.fit_RoBTT                   <- function(object, i){
  
  model              <- object[["models"]][[i]]
  priors             <- model[["priors"]]
  control            <- object[["control"]]
  add_info           <- object[["add_info"]]
  convergence_checks <- object[["convergence_checks"]]
  
  # prepare fit data
  fit_data <- .fit_data(object[["data"]], priors, model[["likelihood"]])
  
  # fit the model
  fit <- .fit_model_and_marglik_RoBTT(
    priors             = priors,
    fit_data           = fit_data,
    likelihood         = model[["likelihood"]],
    control            = control,
    convergence_checks = convergence_checks
  )
  
  model <- c(model, fit)
  
  return(model)
}
.fit_model_and_marglik_RoBTT <- function(priors, fit_data, likelihood, control, convergence_checks){
  
  fit_summary <- NULL
  errors      <- NULL
  warnings    <- NULL
  
  model_call <- list(
    object    = stanmodels[[likelihood]],
    data      = fit_data,
    chains    = control[["chains"]],
    warmup    = control[["warmup"]],
    iter      = control[["iter"]],
    thin      = control[["thin"]],
    cores     = control[["cores"]],
    control   = list(
      adapt_delta   = control[["adapt_delta"]],
      max_treedepth = control[["max_treedepth"]]
    )
  )
  
  if(control[["silent"]]){
    model_call$refresh <- -1
    model_call$open_progress <- FALSE
    model_call$show_messages <- FALSE
  }
  
  if(likelihood == "beta"){
    model_call$init <- lapply(1:control[["chains"]], function(i) {
      list(
        mu      = 0.5,
        sigma2  = 1/12
      )
    })
  }
  
  if(!is.null(control[["seed"]])){
    set.seed(control[["seed"]])
    model_call$seed <- control[["seed"]]
  }
  
  fit <- tryCatch(suppressWarnings(do.call(rstan::sampling, model_call)), error = function(e)e)
  
  # for BayesTools formatting
  attr(fit, "prior_list") <- priors
  
  if(all(class(fit) %in% c("simpleError", "error", "condition"))){
    errors    <- c(errors, fit$message)
    converged <- FALSE
  }
  
  if(!is.null(fit) & !any(class(fit) %in% c("simpleError", "error", "condition"))){
    
    fit_summary <- BayesTools::stan_estimates_table(fit)
    converged   <- TRUE 
    
    if(any(fit_summary[,"ESS"] < convergence_checks[["min_ESS"]])){
      warnings  <- c(warnings, paste0("Minimum effective sample size was low (",  round(min(fit_summary[,"ESS"])), ")."))
    }
    
    if(any(fit_summary[,"R_hat"] > convergence_checks[["max_Rhat"]])){
      warnings  <- c(warnings, paste0("Maximum R-hat was large (",  round(max(fit_summary[,"R_hat"]), 2), ")."))
    }
    
    if(any(rstan::get_divergent_iterations(fit))){
      warnings  <- c(warnings, paste0("There were ", sum(rstan::get_divergent_iterations(fit)), " divergent transitions."))
    }
    rstan::check_divergences(fit)
  }
  
  
  if(!is.null(control$seed)){
    set.seed(control$seed)
  }
  
  marglik <- tryCatch(suppressWarnings(bridgesampling::bridge_sampler(
    samples   = fit,
    maxiter   = control[["bridge_max_iter"]],
    silent    = TRUE)),
    error = function(e)return(e))
  
  # handle errors
  if(any(class(marglik) %in% c("simpleError", "error"))){
    

    errors      <- c(errors, marglik$message)
    marglik     <- .marglik_fail()
    converged   <- FALSE
    
  }else if(is.na(marglik$logml)){
    
    errors      <- c(errors, "not enough iterations")
    marglik     <- .marglik_fail()
    converged   <- FALSE
  }
  
  return(list(
    fit         = fit,
    marglik     = marglik,
    fit_summary = fit_summary,
    errors      = errors,
    warnings    = warnings,
    converged   = converged
  ))
}
.fit_data              <- function(data, priors, likelihood){
  
  data <- .stan_data(data)
  data <- c(data, .stan_distribution("d", priors[["delta"]]))
  data <- c(data, .stan_distribution("r", priors[["rho"]]))
  if(likelihood == "t"){
    data <- c(data, .stan_distribution("nu", priors[["nu"]]))
  }
  
  return(data)
}
.marglik_fail          <- function(){
  marglik        <- NULL
  marglik$logml  <- -Inf
  class(marglik) <- "bridge"
  return(marglik)
}
.fitting_priority      <- function(models){
  
  # model fitting difficulty using the following heuristic:
  # non-normal > heterogeneity > effect
  fitting_difficulty <- sapply(models, function(model){
    
    difficulty <- 0
    
    if(!is.null(model$priors[["delta"]]) && is.prior.simple(model$priors[["delta"]])){
      difficulty <- difficulty + 1
    }
    if(!is.null(model$priors[["rho"]])   && is.prior.simple(model$priors[["rho"]])){
      difficulty <- difficulty + 2
    }
    if(!is.null(model$priors[["nu"]])    && is.prior.simple(model$priors[["nu"]])){
      difficulty <- difficulty + 3
    }
    
    return(difficulty)
  })
  
  return(order(fitting_difficulty, decreasing = TRUE))
}
