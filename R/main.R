#' @title Estimate a Robust Bayesian T-Test
#'
#' @description \code{RoBTT} is used to estimate a Robust Bayesian
#' T-Test. The input either requires the vector of observations for 
#' each group, \code{x1, x2}, or the summary statistics (in case only 
#' the \code{"normal"} likelihood is used).
#'
#' @param x1 vector of observations of the first group
#' @param x2 vector of observations of the second group
#' @param mean1 mean of the first group
#' @param mean2 mean of the first group
#' @param sd1 standard deviation of the first group
#' @param sd2 standard deviation of the first group
#' @param N1 sample size of the first group
#' @param N2 sample size of the first group
#' @param prior_delta prior distributions for the effect size \code{delta} parameter 
#' that will be treated as belonging to the alternative hypothesis. Defaults to \code{
#' prior(distribution = "Cauchy", parameters = list(location = 0, scale = sqrt(2)/2))}.
#' @param prior_rho prior distributions for the precision allocation \code{rho} parameter 
#' that will be treated as belonging to the alternative hypothesis. Defaults to \code{
#' prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))}.
#' @param prior_nu prior distribution for the degrees of freedom + 2 \code{nu}
#' parameter that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{prior(distribution = "exp", parameters = list(rate = 1))}.
#' @param prior_delta_null prior distribution for the \code{delta} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param prior_rho_null prior distribution for the \code{rho} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0.5 (
#' \code{prior(distribution = "point", parameters = list(location = 0.5))}).
#' @param prior_nu_null prior distribution for the \code{nu} parameter
#' that will be treated as belonging to the null hypothesis. Defaults to NULL (
#' \code{NULL}, i.e., normal distribution).
#' @param likelihood types of likelihood to be used. Defaults to \code{c("normal", "t")} 
#' for a normal and t likelihoods.
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
#'
#' @details Add more details
#'
#' Generic [summary.RoBTT()], [print.RoBTT()], and [plot.RoBTT()] functions are
#' provided to facilitate manipulation with the ensemble.
#'
#' @return \code{RoBTT} returns an object of \link[base]{class} \code{"RoBTT"}.
#'
#' @examples \dontrun{
#' # using the example data from XXX
#' }
#'
#' @references
#' \insertAllCited{}
#' @export RoBTT
#' @seealso [summary.RoBTT()], [prior()]
RoBTT <- function(
  x1 = NULL, x2 = NULL,
  mean1 = NULL, mean2 = NULL, sd1 = NULL, sd2 = NULL, N1 = NULL, N2 = NULL,
  
  prior_delta  = prior(distribution = "cauchy",  parameters = list(location = 0, scale = sqrt(2)/2)),
  prior_rho    = prior(distribution = "beta",    parameters = list(alpha = 1, beta = 1)),
  prior_nu     = prior(distribution = "exp",     parameters = list(rate = 1)),
  
  prior_delta_null  = prior(distribution = "spike",  parameters = list(location = 0)),
  prior_rho_null    = prior(distribution = "spike",  parameters = list(location = 0.5)),
  prior_nu_null     = NULL,

  likelihood = c("normal", if(!is.null(prior_nu)) "t"),
  
  chains  = 4, iter = 10000, warmup = 5000, thin = 1, parallel = FALSE,
  control = set_control(), convergence_checks = set_convergence_checks(), 
  
  save = "all", seed = NULL, silent = TRUE, ...){
  
  dots         <- .RoBTT_collect_dots(...)
  object       <- NULL
  object$call  <- match.call()
  object$data  <- .check_data(x1 = x1, x2 = x2, mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, N1 = N1, N2 = N2)

  ### check MCMC settings
  object$control            <- .stan_check_and_list_fit_settings(chains = chains, warmup = warmup, iter = iter, thin = thin,
                                                                 parallel = parallel, cores = chains, silent = silent, seed = seed, control = control)
  object$convergence_checks <- .check_and_list_convergence_checks(convergence_checks)
  
  ### prepare and check the settings
  object$priors      <- .set_priors(prior_delta, prior_rho, prior_nu, prior_delta_null, prior_rho_null, prior_nu_null)
  object$models      <- .get_models(object$priors, likelihood)
  object$add_info$warnings <- c()
  

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
