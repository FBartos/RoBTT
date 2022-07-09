.ensemble_inference  <- function(object){
  
  # modify the default null prior distribution to be spike at Inf
  for(i in seq_along(object[["models"]])){
    if(object[["models"]][[i]]$likelihood == "normal"){
      attr(object[["models"]][[i]][["fit"]], "prior_list")$nu <- prior("spike", parameters = list(location = Inf))
    }
  }
  
  # use only converged models with prior weights > 0 for inference about parameters
  prior_weights <- sapply(object[["models"]], function(model) model[["prior_weights"]])
  models        <- object[["models"]][.get_model_convergence(object) & prior_weights > 0]
  
  # obtain the component type
  effect        <- sapply(models, function(m)!.is_parameter_null(m$priors, "delta"))
  heterogeneity <- sapply(models, function(m)!.is_parameter_null(m$priors, "rho"))
  outliers      <- sapply(models, function(m)!.is_parameter_null(m$priors, "nu"))
  
  # define inference options
  components_present <- c(sum(effect), sum(heterogeneity), sum(outliers)) > 0
  
  components      <- c("Effect", "Heterogeneity", "Outliers")[components_present]
  parameters      <- c("delta", "rho", "nu")[components_present]
  components_null <- list("Effect" = !effect, "Heterogeneity" = !heterogeneity, "Outliers" = !outliers)[components_present]
  parameters_null <- list("mu"     = !effect, "tau"           = !heterogeneity, "nu"       = !outliers)[components_present]
  
  ### get models inference
  inference <- BayesTools::ensemble_inference(
    model_list   = models,
    parameters   = components,
    is_null_list = components_null,
    conditional  = FALSE
  )
  # deal with the possibility of only null models models
  if(all(sapply(components_null, all))){
    inference_conditional <- NULL
  }else{
    inference_conditional <- BayesTools::ensemble_inference(
      model_list   = models,
      parameters   = components[!sapply(components_null, all)],
      is_null_list = components_null[!sapply(components_null, all)],
      conditional  = TRUE
    )
  }
  
  
  ### get model-averaged posteriors
  posteriors <- BayesTools::mix_posteriors(
    model_list   = models,
    parameters   = parameters,
    is_null_list = parameters_null,
    seed         = object$add_info[["seed"]],
    conditional  = FALSE
  )
  # deal with the possibility of only null models models
  if(all(sapply(components_null, all))){
    posteriors_conditional <- NULL
  }else{
    posteriors_conditional <- BayesTools::mix_posteriors(
      model_list   = models,
      parameters   = parameters[!sapply(parameters_null, all)],
      is_null_list = parameters_null[!sapply(parameters_null, all)],
      seed         = object$add_info[["seed"]],
      conditional  = TRUE
    )
  }
  
  ### get posteriors for the mu and sigma
  parameters_est      <- c("mu_i", "sigma_i")
  parameters_est_null <- list("mu_i" = rep(FALSE, length(models)), "sigma_i" = rep(FALSE, length(models)))
  
  # add fake vector priors 
  for(i in seq_along(models)){
    attr(models[[i]][["fit"]], "prior_list") <- c(
      attr(models[[i]][["fit"]], "prior_list"),
      list("mu_i"    = BayesTools::prior("mnormal", list(mean = 0, sd = 1, K = 2))),
      list("sigma_i" = BayesTools::prior("mnormal", list(mean = 0, sd = 1, K = 2)))
    )
  }

  posteriors_est <- BayesTools::mix_posteriors(
    model_list   = models,
    parameters   = parameters_est,
    is_null_list = parameters_est_null,
    seed         = object$add_info[["seed"]],
    conditional  = FALSE
  )
  
  # remove "_i" from the name
  colnames(posteriors_est[["mu_i"]])    <- paste0("mu[", 1:2, "]")
  colnames(posteriors_est[["sigma_i"]]) <- paste0("sigma[", 1:2, "]")

  
  # return the results
  output <- list(
    inference              = inference,
    inference_conditional  = inference_conditional,
    posteriors             = posteriors,
    posteriors_conditional = posteriors_conditional,
    posteriors_est         = posteriors_est
  )
  
  return(output)
}
.compute_coeficients <- function(RoBTT){
  return(c(
    "delta" = if(length(RoBTT$posteriors[["delta"]]) != 0)  mean(RoBTT$posteriors[["delta"]]),
    "rho"   = if(length(RoBTT$posteriors[["rho"]])   != 0)  mean(RoBTT$posteriors[["rho"]])
  ))
}



