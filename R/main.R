#' @title Estimate a Robust Bayesian T-Test
#'
#' @description \code{RoBTT} is used to estimate a robust Bayesian
#' t-test or truncated Bayesian t-test (if \code{truncation} is used). 
#' The input either requires the vector of observations for 
#' each group, \code{x1, x2}, or the summary statistics (only if the normal 
#' likelihood models are used).
#'
#' @param x1 vector of observations of the first group
#' @param x2 vector of observations of the second group
#' @param mean1 mean of the first group
#' @param mean2 mean of the first group
#' @param sd1 standard deviation of the first group
#' @param sd2 standard deviation of the first group
#' @param N1 sample size of the first group
#' @param N2 sample size of the first group
#' @param truncation an optional list specifying truncation applied to the data. 
#' Defaults to \code{NULL}, i.e., no truncation was applied and the full likelihood is 
#' applied. Alternative the truncation can be specified via a named list with:
#' \describe{
#'   \item{\code{"x"}}{where \code{x} is a vector of two values specifying the lower 
#'   and upper truncation points common across the groups}
#'   \item{\code{"x1"} and \code{"x2"}}{where \code{x1} is a vector of two values specifying 
#'   the lower and upper truncation points for the first group and \code{x2} is a vector of
#'   two values specifying the lower and upper truncation points for the second group.}
#'   \item{\code{"sigma"}}{where \code{sigma} corresponds to the number of standard deviations
#'   from the common mean where the truncation points should be set.}
#'   \item{\code{"sigma1"} and \code{"sigma2"}}{where \code{sigma1} corresponds to the number of
#'   standard deviations from the mean of the first group where the truncation points should be set
#'   and \code{sigma2} corresponds to the number of standard deviations from the mean of the second
#'   group where the truncation points should be set.}
#' }
#' @param prior_delta prior distributions for the effect size \code{delta} parameter 
#' that will be treated as belonging to the alternative hypothesis. Defaults to \code{
#' prior(distribution = "Cauchy", parameters = list(location = 0, scale = sqrt(2)/2))}.
#' @param prior_rho prior distributions for the precision allocation \code{rho} parameter 
#' that will be treated as belonging to the alternative hypothesis. Defaults to \code{
#' prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))}.
#' @param prior_nu prior distribution for the degrees of freedom + 2 \code{nu}
#' parameter that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{prior(distribution = "exp", parameters = list(rate = 1))} if no 
#' \code{truncation} is specified. If \code{truncation} is specified, the default is
#' \code{NULL} (i.e., use only normal likelihood).
#' @param prior_delta_null prior distribution for the \code{delta} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param prior_rho_null prior distribution for the \code{rho} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0.5 (
#' \code{prior(distribution = "point", parameters = list(location = 0.5))}).
#' @param prior_nu_null prior distribution for the \code{nu} parameter
#' that will be treated as belonging to the null hypothesis. Defaults to \code{prior_none} (
#' (i.e., normal likelihood)).
#' @param prior_mu prior distribution for the grand mean parameter. Defaults to \code{NULL} 
#' which sets Jeffreys prior for the grand mean in case of no truncation or an unit Cauchy 
#' prior distributions for the grand mean in case of truncation (which greatly improves 
#' sampling efficiency).
#' @param prior_sigma prior distribution for the grand variance parameter. Defaults to \code{NULL}
#' which sets Jeffreys prior for the standard deviation in case of no truncation or an exponential prior
#' distribution for the standard deviation in case of truncation (which greatly improves sampling efficiency).
#' @param chains a number of chains of the MCMC algorithm.
#' @param iter a number of sampling iterations of the MCMC algorithm.
#' Defaults to \code{10000}, with a minimum of \code{4000}.
#' @param warmup a number of warmup  iterations of the MCMC algorithm.
#' Defaults to \code{5000}.
#' @param thin a thinning of the chains of the MCMC algorithm. Defaults to
#' \code{1}.
#' @param parallel whether the individual models should be fitted in parallel.
#' Defaults to \code{FALSE}. The implementation is not completely stable
#' and might cause a connection error.
#' @param control allows to pass control settings with the
#' [set_control()] function. See \code{?set_control} for
#' options and default settings.
#' @param convergence_checks automatic convergence checks to assess the fitted
#' models, passed with [set_convergence_checks()] function. See
#' \code{?set_convergence_checks} for options and default settings.
##' @param save whether all models posterior distributions should be kept
#' after obtaining a model-averaged result. Defaults to \code{"all"} which
#' does not remove anything. Set to \code{"min"} to significantly reduce
#' the size of final object, however, some model diagnostics and further
#' manipulation with the object will not be possible.
#' @param seed a seed to be set before model fitting, marginal likelihood
#' computation, and posterior mixing for reproducibility of results. Defaults
#' to \code{NULL} - no seed is set.
#' @param silent whether all print messages regarding the fitting process
#' should be suppressed. Defaults to \code{TRUE}. Note that \code{parallel = TRUE}
#' also suppresses all messages.
#' @param ... additional arguments.
#'
#' @details 
#' See \insertCite{maier2022bayesian;textual}{RoBTT} for more details 
#' regarding the robust Bayesian t-test methodology and the corresponding 
#' vignette (\href{../doc/Introduction_to_RoBTT.html}{\code{vignette("Introduction_to_RoBTT", package = "RoBTT")}}).
#' 
#' See \insertCite{godmann2024how;textual}{RoBTT} for more details 
#' regarding the truncated Bayesian t-test methodology and the corresponding 
#' vignette (\href{../doc/Truncated_t_test.html}{\code{vignette("Truncated_t_test", package = "RoBTT")}}).
#'
#' Generic [summary.RoBTT()], [print.RoBTT()], and [plot.RoBTT()] 
#' functions are provided to facilitate manipulation with the ensemble.
#' 
#'
#' @return \code{RoBTT} returns an object of \link[base]{class} \code{"RoBTT"}.
#'
#' @examples \dontrun{
#' # using the example data from Darwin
#' data("fertilization", package = "RoBTT")
#' fit <- RoBTT(
#'   x1       = fertilization$Self,
#'   x2       = fertilization$Crossed,
#'   prior_delta = prior("cauchy", list(0, 1/sqrt(2))),
#'   prior_rho   = prior("beta",   list(3, 3)),
#'   seed        = 1, 
#'   chains      = 1,
#'   warmup      = 1000,
#'   iter        = 2000,
#'   control     = set_control(adapt_delta = 0.95)
#' )
#'
#' # summary can provide many details about the model
#' summary(fit)
#' }
#' 
#' @references
#' \insertAllCited{}
#' @export RoBTT
#' @seealso [summary.RoBTT()], [prior()]
RoBTT <- function(
  x1 = NULL, x2 = NULL,
  mean1 = NULL, mean2 = NULL, sd1 = NULL, sd2 = NULL, N1 = NULL, N2 = NULL,
  truncation = NULL,
  
  prior_delta  = prior(distribution = "cauchy",  parameters = list(location = 0, scale = sqrt(2)/2)),
  prior_rho    = prior(distribution = "beta",    parameters = list(alpha = 1, beta = 1)),
  prior_nu     = if(is.null(truncation)) prior(distribution = "exp", parameters = list(rate = 1)),
  
  prior_delta_null  = prior(distribution = "spike",  parameters = list(location = 0)),
  prior_rho_null    = prior(distribution = "spike",  parameters = list(location = 0.5)),
  prior_nu_null     = prior_none(),
  
  prior_mu = NULL, prior_sigma = NULL,
  
  chains  = 4, iter = 10000, warmup = 5000, thin = 1, parallel = FALSE,
  control = set_control(), convergence_checks = set_convergence_checks(), 
  
  save = "all", seed = NULL, silent = TRUE, ...){
  
  dots         <- .RoBTT_collect_dots(...)
  object       <- NULL
  object       <- NULL
  object$call  <- match.call()
  object$data  <- .check_data(x1 = x1, x2 = x2, mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, N1 = N1, N2 = N2, truncation = truncation)

  ### check MCMC settings
  object$control            <- .stan_check_and_list_fit_settings(chains = chains, warmup = warmup, iter = iter, thin = thin,
                                                                 parallel = parallel, cores = chains, silent = silent, seed = seed, control = control)
  object$convergence_checks <- .check_and_list_convergence_checks(convergence_checks)
  
  ### prepare and check the settings
  object$priors      <- .set_priors(prior_delta, prior_rho, prior_nu, prior_delta_null, prior_rho_null, prior_nu_null, prior_mu, prior_sigma, object$data, !is.null(truncation))
  object$models      <- .get_models(object$priors)
  object$add_info    <- list(
    warnings         = NULL,
    seed             = seed,
    save             = save
  )
  

  ### fit the models and compute marginal likelihoods
  if(!object$control[["parallel"]]){
    
    if(dots[["is_JASP"]]){
      .JASP_progress_bar_start(length(object[["models"]]))
    }
    
    for(i in seq_along(object[["models"]])){
      object$models[[i]] <- .fit_RoBTT(object, i)
      if(dots[["is_JASP"]]){
        .JASP_progress_bar_tick()
      }
    }
    
  }else{
    
    fitting_order <- .fitting_priority(object[["models"]])
    
    cl <- parallel::makePSOCKcluster(floor(RoBTT.get_option("max_cores") / object$control[["chains"]]))
    parallel::clusterEvalQ(cl, {library("RoBTT")})
    parallel::clusterExport(cl, "object", envir = environment())
    object$models <- parallel::parLapplyLB(cl, fitting_order, .fit_RoBTT, object = object)[order(fitting_order)]
    parallel::stopCluster(cl)
    
  }

  # deal with non-converged the converged models
  object$add_info$converged <- .get_model_convergence(object)

  # create ensemble only if all models converged
  if(all(object$add_info$converged)){
    ### compute the model-space results
    object$models        <- BayesTools::models_inference(object[["models"]])
    object$RoBTT         <- .ensemble_inference(object)
    object$coefficients  <- .compute_coeficients(object[["RoBTT"]])
  }else{
    stop(paste0("The following models failed the converge: ", paste(which(!object$add_info$converged), collapse = ", "), "."))
  }


  ### collect and print errors and warnings
  object$add_info[["errors"]]   <- c(object$add_info[["errors"]],   .get_model_errors(object))
  object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], .get_model_warnings(object))
  .print_errors_and_warnings(object)
  
  
  ### remove model posteriors if asked to
  if(save == "min"){
    object <- .remove_model_posteriors(object)
    object <- .remove_model_margliks(object)
  }

  class(object) <- "RoBTT"
  return(object)
}


#' @title Updates a fitted RoBTT object
#'
#' @description \code{update.RoBTT} can be used to
#' \enumerate{
#'   \item{change the prior odds of fitted models by specifying a vector
#'   \code{prior_weights} of the same length as the fitted models,}
#'   \item{refitting models that failed to converge with updated settings
#'   of control parameters,}
#'   \item{or changing the convergence criteria and recalculating the ensemble
#'   results by specifying new \code{control} argument and setting
#'   \code{refit_failed == FALSE}.}
#' }
#'
#' @param object a fitted RoBTT object
#' @param prior_weights either a single value specifying prior model weight
#' of a newly specified model using priors argument, or a vector of the
#' same length as already fitted models to update their prior weights.
#' @param refit_failed whether failed models should be refitted. Relevant only
#' \code{prior_weights} are not supplied. Defaults to \code{TRUE}.
#' @inheritParams RoBTT
#' @param ... additional arguments.
#'
#' @details See [RoBTT()] for more details.
#'
#' @return \code{RoBTT} returns an object of class 'RoBTT'.
#'
#' @seealso [RoBTT()], [summary.RoBTT()], [prior()], [check_setup()]
#' @export
update.RoBTT <- function(object, refit_failed = TRUE, prior_weights = NULL,
                         chains  = NULL, iter = NULL, warmup = NULL, thin = NULL, parallel = NULL,
                         control = NULL, convergence_checks = NULL,
                         save = "all", seed = NULL, silent = TRUE, ...){
  
  BayesTools::check_bool(refit_failed, "refit_failed")
  dots <- .RoBTT_collect_dots(...)

  if(object$add_info$save == "min")
    stop("Models cannot be updated because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' in while fitting the model (see ?RoBMA for more details).")
  
  if(!is.null(prior_weights)){
    
    what_to_do <- "update_prior_weights"
    if(length(prior_weights) != length(object$models))
      stop("The number of newly specified prior odds does not match the number of models. See '?update.RoBTT' for more details.")
    for(i in 1:length(object$models)){
      object$models[[i]]$prior_weights     <- prior_weights[i]
      object$models[[i]]$prior_weights_set <- prior_weights[i]
    }
    
  }else if(refit_failed & any(!.get_model_convergence(object, include_warning = TRUE))){
    
    what_to_do <- "refit_failed_models"
    
  }else{
    
    what_to_do <- "update_settings"
    
  }
  
  
  ### update control settings if any change is specified
  object[["control"]]            <- .update_fit_control(object[["control"]], chains = chains, warmup = warmup, iter = iter, thin = thin, parallel = parallel, cores = chains, silent = silent, seed = seed, control = control)
  object[["convergence_checks"]] <- .update_convergence_checks(object[["convergence_checks"]], convergence_checks)
  
  
  ### clean errors and warnings
  object$add_info[["errors"]]   <- NULL
  object$add_info[["warnings"]] <- NULL
  
  
  ### do the stuff
  if(what_to_do == "refit_failed_models"){
    
    models_to_update <- seq_along(object$models)[!.get_model_convergence(object, include_warning = TRUE)]
    
    if(!object$control[["parallel"]]){
      
      if(dots[["is_JASP"]]){
        .JASP_progress_bar_start(length(models_to_update))
      }
      
      for(i in models_to_update){
        object$models[[i]] <- .fit_RoBTT(object, i)
        if(dots[["is_JASP"]]){
          .JASP_progress_bar_tick()
        }
      }
      
    }else{
      
      cl <- parallel::makePSOCKcluster(floor(RoBTT.get_option("max_cores") / object$control[["chains"]]))
      parallel::clusterEvalQ(cl, {library("RoBTT")})
      parallel::clusterExport(cl, "object", envir = environment())
      object$models[models_to_update] <- parallel::parLapplyLB(cl, models_to_update, .fit_RoBTT, object = object)
      parallel::stopCluster(cl)
      
    }
    
  }
  
  
  # create ensemble only if at least one model converged
  if(all(object$add_info$converged)){
    ### compute the model-space results
    object$models        <- BayesTools::models_inference(object[["models"]])
    object$RoBTT         <- .ensemble_inference(object)
    object$coefficients  <- .compute_coeficients(object[["RoBTT"]])
  }else{
    stop(paste0("The following models failed the converge: ", paste(which(!object$add_info$converged), collapse = ", "), "."))
  }
  
  
  ### collect and print errors and warnings
  object$add_info[["errors"]]   <- c(object$add_info[["errors"]],   .get_model_errors(object))
  object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], .get_model_warnings(object))
  .print_errors_and_warnings(object)
  
  
  ### remove model posteriors if asked to
  if(save == "min"){
    object <- .remove_model_posteriors(object)
    object <- .remove_model_margliks(object)
  }
  
  return(object)
}

