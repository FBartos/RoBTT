.remove_model_posteriors   <- function(object){
  for(i in seq_along(object[["models"]])){
    if(inherits(object$models[[i]][["fit"]], "stanfit")){
      object$models[[i]]$fit <- NULL
    }
  }
  return(object)
}
.remove_model_margliks     <- function(object){
  for(i in seq_along(object[["models"]])){
    if(inherits(object$models[[i]][["marglik"]], "bridge")){
      object$models[[i]]$marglik[["q11"]] <- NULL
      object$models[[i]]$marglik[["q12"]] <- NULL
      object$models[[i]]$marglik[["q21"]] <- NULL
      object$models[[i]]$marglik[["q22"]] <- NULL
    }
  }
  return(object)
}
.print_errors_and_warnings <- function(object, max_print = 5){
  
  errors_and_warnings <- .collect_errors_and_warnings(object, max_print = max_print)
  
  for(i in seq_along(errors_and_warnings))
    warning(errors_and_warnings[i], immediate. = TRUE, call. = FALSE)
  
  return(invisible(errors_and_warnings))
}
.shorten_warnings    <- function(warnings, n_warnings = 5){
  if(is.null(warnings)){
    return(NULL)
  }else if(length(warnings) <= n_warnings){
    return(warnings)
  }else{
    return(c(warnings[1:n_warnings], paste0("There were another ", length(warnings) - n_warnings - 1, " warnings. To see all warnings call 'check_RoBMA(fit)'.")))
  }
}
.shorten_errors      <- function(errors, n_errors = 5){
  if(is.null(errors)){
    return(NULL)
  }else if(length(errors) <= n_errors){
    return(errors)
  }else{
    return(c(errors[1:n_errors], paste0("There were another ", length(errors) - n_errors - 1, " errors. To see all errors call 'check_RoBMA(fit)'.")))
  }
}
.convergence_warning <- function(object){
  if(any(!.get_model_convergence(object))){
    return(paste0(sum(!.get_model_convergence(object)), ifelse(sum(!.get_model_convergence(object)) == 1, " model", " models"), " failed to converge."))
  }else{
    return(NULL)
  }
}
.collect_errors_and_warnings <- function(object, max_print = 5){
  
  short_warnings <- .shorten_warnings(object$add_info[["warnings"]], max_print)
  short_errors   <- .shorten_errors(object$add_info[["errors"]],     max_print)
  conv_warning   <- .convergence_warning(object)
  
  return(c(short_warnings, short_errors, conv_warning))
}
.get_model_convergence       <- function(object){
  return(sapply(object[["models"]], function(model) if(is.null(model[["converged"]])) FALSE else model[["converged"]]))
}
.get_model_warnings          <- function(object){
  return(unlist(sapply(seq_along(object[["models"]]), function(i){
    if(!is.null(object[["models"]][[i]][["warnings"]])){
      paste0("Model (", i, "): ", object[["models"]][[i]][["warnings"]])
    }
  })))
}
.get_model_errors            <- function(object){
  return(unlist(sapply(seq_along(object[["models"]]), function(i){
    if(!is.null(object[["models"]][[i]][["errors"]])){
      paste0("Model (", i, "): ", object[["models"]][[i]][["errors"]])
    }
  })))
}
.is_parameter_null <- function(priors, par){
  return(if(is.null(priors[[par]])) TRUE else priors[[par]][["is_null"]])
}
.is_model_constant <- function(priors){
  
  constant <- NULL
  for(par in c("mu", "tau", "omega", "sigma")){
    if(!is.null(priors[[par]])){
      constant <- c(constant, is.prior.point(priors[[par]]))
    }
  }
  
  constant <- all(constant)
  
  return(constant)
}
.get_no_support    <- function(models, par){
  
  no_support  <- NULL
  
  all_support <- sapply(models, function(m) m$priors[[par]]$truncation, simplify = FALSE)
  all_support <- do.call(rbind.data.frame, all_support)
  
  if(!is.null(all_support)){
    
    # start
    if(!is.infinite(min(all_support$lower))){
      no_support <- c(no_support, list(list(lower = -Inf, upper = min(all_support$lower))))
      temp_end   <- min(all_support$lower)
    }else{
      temp_end   <- -Inf
    }
    
    # the middle
    all_support <- all_support[order(all_support$lower),]
    for(i in 1:nrow(all_support)){
      
      # prolong the current coverage
      if(all_support$lower[i] <= temp_end & all_support$upper[i] > temp_end){
        temp_end <- all_support$upper[i]
        next
      }
      
      # detect the gap
      if(all_support$lower[i] > temp_end){
        no_support <- c(no_support, list(list(lower = temp_end, upper = all_support$lower[i])))
        temp_end   <- all_support$lower[i]
      }
      
    }
    
    # the upper part
    if(!is.infinite(max(all_support$upper)))no_support <- c(no_support, list(list(lower = max(all_support$upper), upper = Inf)))
  }
  
  return(no_support)
}
.stan.summary      <- function(fit){
  
  summary_fit   <- rstan::summary(fit)$summary
  
  return(summary_fit)
}
