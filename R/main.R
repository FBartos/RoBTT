#' @title Estimate a Robust Bayesian T-Test
#'
#' @description \code{RoBTT} is used to estimate a Robust Bayesian
#' T-Test. 
#'
#' @param priors_mu list of prior distributions for the \code{mu} parameter that
#' will be treated as belonging to the alternative hypothesis. Defaults to \code{
#' prior(distribution = "normal",   parameters = list(mean = 0, sd = 1))}.
#' @param priors_tau list of prior distributions for the \code{tau} parameter that
#' will be treated as belonging to the alternative hypothesis. Defaults to \code{
#' prior(distribution = "invgamma", parameters = list(shape = 1, scale = .15))}.
#' @param priors_omega list of prior weight functions for the \code{omega}
#' parameter that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{list(
#' prior(distribution = "two.sided", parameters = list(alpha = c(1, 1),     steps = c(.05)),      prior_odds = 1/2),
#' prior(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),  steps = c(.05, .10)), prior_odds = 1/2)
#' )}.
#' @param priors_mu_null list of prior distributions for the \code{mu} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param priors_tau_null list of prior distributions for the \code{tau} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param priors_omega_null list of prior weight functions for the \code{omega} parameter
#' that will be treated as belonging to the null hypothesis. Defaults to point
#' distribution with location at 1 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param chains a number of chains of the MCMC algorithm.
#' @param iter a number of sampling iterations of the MCMC algorithm.
#' Defaults to \code{10000}, with a minimum of \code{4000}.
#' @param warmup a number of warmup  iterations of the MCMC algorithm.
#' Defaults to \code{5000}.
#' @param thin a thinning of the chains of the MCMC algorithm. Defaults to
#' \code{1}.
#' @param control a list of additional arguments for the MCMC algorithm.
#' Possible options are:
#' \describe{
#'   \item{adapt_delta}{Defaults to \code{.80}.}
#'   \item{max_treedepth}{Defaults to \code{15}.}
#'   \item{bridge_max_iter}{Maximum number of iterations for the
#'   \link[bridgesampling]{bridge_sampler} function. Defaults to \code{10000}}
#'   \item{allow_max_rhat}{Maximum allowed Rhat for a model to be taken into
#'   consideration. Model will be removed from the ensemble if it fails to
#'   achieve the set Rhat. Defaults to \code{NULL} - no model will be removed
#'   based on Rhat.}
#'   \item{allow_min_ESS}{Minimum allowed ESS for a model to be taken into
#'   consideration. Model will be removed from the ensemble if it fails to
#'   achieve the set ESS. Defaults to \code{NULL} - no model will be removed
#'   based on ESS.}
#'   \item{silent}{Whether all fitting messages should be suppressed. Defaults
#'   to \code{FALSE}. Ideal for getting rid of the "full precision may not have
#'   been achieved in pnt{final}'" warning that cannot be suppressed in any
#'   other way.}
#'   \item{cores}{Maximum number of cores to be used for parallel computation. If
#'   \code{parallel == TRUE}, the default number is equal to number of cores - 1,
#'   and 1 (no parallel processing otherwise).}
#' }
#' @param seed a seed to be set before model fitting, marginal likelihood
#' computation, and posterior mixing for exact results reproducibility. Defaults
#' to \code{NULL} - no seed is set.
#' @param parallel whether the individual models should be fitted in parallel.
#' Defaults to \code{FALSE}. The \code{cores} argument within the \code{control}
#' list will overwrite the setting if specified to a number higher than 1.
#'
#' @details Add more details
#'
#' Generic [summary.RoBTT()], [print.RoBTT()], and [plot.RoBTT()] functions are
#' provided to facilitate manipulation with the ensemble. A visual check of the
#' individual model diagnostics can be obtained using the [diagnostics()] function.
#' The fitted model can be further updated or modified by [update.RoBTT()] function.
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
#' @seealso [summary.RoBTT()], [update.RoBTT()], [prior()], [check_setup()]
RoBTT <- function(
  x1 = NULL, x2 = NULL,
  mean1 = NULL, mean2 = NULL, sd1 = NULL, sd2 = NULL, N1 = NULL, N2 = NULL,
  
  prior_d   = prior(distribution = "cauchy",  parameters = list(location = 0, scale = sqrt(2)/2)),
  prior_r   = prior(distribution = "beta",    parameters = list(alpha = 1, beta = 1)),
  prior_nu  = prior(distribution = "exp",     parameters = list(rate = 1)),
  prior_p = prior(distribution = "normal",  parameters = list(mean = 0, sd = 1)),
  prior_d_null   = prior(distribution = "spike",  parameters = list(location = 0)),
  prior_r_null   = prior(distribution = "spike",  parameters = list(location = 0.5)),
  prior_nu_null  = NULL,
  prior_p_null   = prior(distribution = "spike",  parameters = list(location = 0)),
  likelihood = c("normal", if(!is.null(prior_nu)) "t"),
  chains  = 4, iter = 10000, warmup = 5000, thin = 1, parallel = FALSE,
  control = NULL, seed = NULL){
  
  object       <- NULL
  object$call  <- match.call()
  object$data  <- list(x1 = x1, x2 = x2, mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, N1 = N1, N2 = N2)

  ### prepare and check the settings
  object$priors   <- .set_priors(prior_d, prior_r, prior_nu, prior_p, prior_d_null, prior_r_null, prior_nu_null, prior_p_null)
  object$models   <- .get_models(object$priors, likelihood)
  object$control  <- .set_control(control, chains, iter, warmup, thin, seed, parallel)
  object$add_info$warnings <- c()


  ### fit the models and compute marginal likelihoods
  if(object$control$cores < 2*object$control$chains){

    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    for(i in 1:length(object$models)){
      object$models[[i]] <- .fit_RoBTT_wrap(object, i)
    }

  }else{

    cl <- parallel::makePSOCKcluster(floor(object$control$cores / object$control$chains))
    parallel::clusterEvalQ(cl, {library("RoBTT")})
    parallel::clusterExport(cl, "object", envir = environment())
    object$models <- parallel::clusterApplyLB(cl, 1:length(object$models), .fit_RoBTT_wrap, object = object)
    parallel::stopCluster(cl)

  }
  
  # deal with non-converged the converged models
  object$add_info$converged <- .get_converged_models(object)


  # create ensemble only if at least one model converges
  if(any(object$add_info$converged)){
    ### compute the model-space results
    object$RoBTT         <- .model_inference(object)
    object$coefficients  <- .compute_coeficients(object$RoBTT)
  }


  ### add warnings
  # object$add_info$warnings <- c(object$add_info$warnings, .model_refit_warnings(sapply(1:length(object$models), function(i)object$models[[i]]$metadata, simplify = FALSE)))
  # object$add_info$warnings <- c(object$add_info$warnings, .model_convergence_warnings(object))


  ### remove model posteriors if asked to
  # if(save == "min"){
  #   for(i in 1:length(object$models)){
  #     if(length(object$models[[1]]$fit) != 1){
  #       object$models[[i]]$fit$mcmc <- NULL
  #     }
  #   }
  # }


  ### print warnings
  if(!is.null(object$add_info$warnings)){
    for(w in object$add_info$warnings)warning(w)
  }
  if(sum(!object$add_info$converged) > 0)warning(paste0(sum(!object$add_info$converged), ifelse(sum(!object$add_info$converged) == 1, " model", " models"), " failed to converge."))

  class(object) <- "RoBTT"
  return(object)
}

### fitting function
.fit_RoBTT             <- function(object, i){

  model      <- object$models[[i]]
  priors     <- model$priors
  control    <- object$control
  add_info   <- object$add_info
  refit_info <- NULL


  # prepare fit data
  fit_data <- .fit_data(object$data, priors, model$likelihood)
  
  # fit the model
  fit_out  <- .fit_model_RoBTT_wrap(fit_data, model$likelihood, control)
  
  model <- c(model, fit_out)
  model$metadata$i <- i

  return(model)
}
.fit_model_RoBTT_wrap  <- function(fit_data, likelihood, control){
  if(control$silent){
    fit <- callr::r(
      .fit_model_RoBTT,
      args = list(
        fit_data    = fit_data,
        likelihood  = likelihood,
        control     = control
      )
    )
  }else{
    fit <- .fit_model_RoBTT(
      fit_data    = fit_data,
      likelihood  = likelihood,
      control     = control
    )
  }
  return(fit)
}
.fit_model_RoBTT       <- function(fit_data, likelihood, control){
  requireNamespace("RoBTT")
  
  metadata    <- list()
  fit_summary <- NULL
  
  model_call <- list(
    object          = stanmodels[[likelihood]],
    data            = fit_data,
    chains          = control$chains,
    warmup          = control$warmup,
    iter            = control$iter,
    thin            = control$thin,
    cores           = control$cores
  )
  
  if(likelihood %in% c("beta", "beta01")){
    model_call$init <- lapply(1:control$chains, function(i) {
      list(
        mu      = 0.5,
        sigma2  = 1/12
      )
    })
  }
  
  if(!is.null(control$seed)){
    set.seed(control$seed)
    model_call$seed <- control$seed
  }
  
  fit <- tryCatch(do.call(rstan::sampling, model_call), error = function(e)e)
  
  if(all(class(fit) %in% c("simpleError", "error", "condition"))){
    metadata$refit_info <- fit$message
  }
  
  if(!is.null(fit) & !any(class(fit) %in% c("simpleError", "error", "condition"))){
    fit_summary <- .stan.summary(fit)
  }
  
  
  if(!is.null(control$seed))set.seed(control$seed)
  marg_lik <- tryCatch(suppressWarnings(bridgesampling::bridge_sampler(
    samples   = fit,
    maxiter   = control$bridge_max_iter,
    silent    = TRUE)),
    error = function(e)return(e))
  
  # handle errors
  if(any(class(marg_lik) %in% c("simpleError", "error"))){
    
    metadata$marg_lik <- marg_lik$message
    marg_lik <- .marglik_fail()
    
  }else if(is.na(marg_lik$logml)){
    
    metadata$marg_lik <- "not enough iterations"
    marg_lik <- .marglik_fail()
    
  }

  return(list(
    fit         = fit,
    marg_lik    = marg_lik,
    fit_summary = fit_summary,
    metadata    = metadata
  ))
}
.marglik_RoBTT         <- function(object, i){

  model    <- object$models[[i]]
  fit      <- model$fit
  control  <- object$control

  if(!is.null(control$seed))set.seed(control$seed)
  marg_lik        <- tryCatch(suppressWarnings(bridgesampling::bridge_sampler(
    samples = fit,
    maxiter = control$bridge_max_iter,
    silent  = TRUE)),
    error = function(e)return(e))

  # handle errors
  if(any(class(marg_lik) %in% c("simpleError", "error"))){

    model$metadata$marg_lik <- marg_lik$message
    marg_lik <- .marglik_fail()

  }else if(is.na(marg_lik$logml)){

    model$metadata$marg_lik <- "not enough iterations"
    marg_lik <- .marglik_fail()

  }

  model$marg_lik <- marg_lik

  return(model)
}
.fit_data              <- function(data, priors, likelihood){

  data <- .stan_data(data$x1, data$x2, data$mean1, data$mean2, data$sd1, data$sd2, data$N1, data$N2)
  data <- c(data, .stan_distribution("d", priors[["d"]]))
  data <- c(data, .stan_distribution("r", priors[["r"]]))
  if(likelihood == "t"){
    data <- c(data, .stan_distribution("nu", priors[["nu"]]))
  }
  if(likelihood %in% c("gamma0", "lognormal0")){
    data <- c(data, .stan_distribution("p", priors[["p"]]))
  }
  
  return(data)
}
.marglik_fail          <- function(){
  marg_lik        <- NULL
  marg_lik$logml  <- -Inf
  class(marg_lik) <- "bridge"
  return(marg_lik)
}
.fit_RoBTT_wrap        <- function(object, i){

  object$models[[i]] <- .fit_RoBTT(object, i)
  object$models[[i]] <- .marglik_RoBTT(object, i)
  if(!is.null(object$control$progress_tick))eval(parse(text = object$control$progress_tick))

  return(object$models[[i]])
}

### model inference functions
.model_inference            <- function(object, n_samples = 10000){

  models      <- object$models
  data        <- object$data
  add_info    <- object$add_info
  converged   <- object$add_info$converged
  seed        <- object$control$seed
  likelihoods <- sapply(models, function(m) m$likelihood)
  
  # extract marginal likelihoods
  marg_liks <- sapply(models, function(x)x$marg_lik$logml)

  # determine the type of the models
  mm_d  <- sapply(models, function(m)!.is_parameter_null(m$priors, "d"))
  mm_r  <- sapply(models, function(m)!.is_parameter_null(m$priors, "r"))
  mm_nu <- sapply(models, function(m)!.is_parameter_null(m$priors, "nu"))
  mm_p  <- sapply(models, function(m)!.is_parameter_null(m$priors, "p"))

  parameters <- c("d", "r", "nu", "p")[c(sum(mm_d), sum(mm_r), sum(mm_nu), sum(mm_p)) > 0]
  
  # extract model weights
  prior_weights_all  <- sapply(models, function(m)m$prior_odds)
  prior_weights_d    <- ifelse(mm_d,  prior_weights_all, 0)
  prior_weights_r    <- ifelse(mm_r,  prior_weights_all, 0)
  prior_weights_nu   <- ifelse(mm_nu, prior_weights_all, 0)
  prior_weights_p    <- ifelse(mm_p,  prior_weights_all, 0)
  
  # standardize model weights
  prior_weights_all  <- prior_weights_all  / sum(prior_weights_all)
  prior_weights_d    <- prior_weights_d    / sum(prior_weights_d)
  prior_weights_r    <- prior_weights_r    / sum(prior_weights_r)
  prior_weights_nu   <- prior_weights_nu   / sum(prior_weights_nu)
  prior_weights_p    <- prior_weights_p    / sum(prior_weights_p)


  ### compute model weights
  # overall
  weights_all   <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_all)
  if(any(mm_d) & all(!is.nan(prior_weights_d))){
    weights_d  <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_d)
  }else{
    weights_d <- NULL
  }
  if(any(mm_r) & all(!is.nan(prior_weights_r))){
    weights_r  <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_r)
  }else{
    weights_r <- NULL
  }
  if(any(mm_nu) & all(!is.nan(prior_weights_nu))){
    weights_nu  <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_nu)
  }else{
    weights_nu <- NULL
  }
  if(any(mm_p) & all(!is.nan(prior_weights_p))){
    weights_p  <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_p)
  }else{
    weights_p <- NULL
  }
  weights_list <- list(
    d  = weights_d,
    r  = weights_r,
    nu = weights_nu,
    p  = weights_p
  )
  
  
  ### compute inclusion BFs
  BF_effect        <- .inclusion_BF(prior_weights_all, weights_all, mm_d)
  BF_heterogeneity <- .inclusion_BF(prior_weights_all, weights_all, mm_r)
  BF_outliers      <- .inclusion_BF(prior_weights_all, weights_all, mm_nu)
  BF_proportion    <- .inclusion_BF(prior_weights_all, weights_all, mm_p)

  
  ### sample and mix the individual posteriors
  if(!is.null(seed))set.seed(seed)
  samples <- list()
  for(par in parameters){
    samples$averaged[[par]]    <- .mix_samples(models, weights_all, converged, par, n_samples, seed)
  }
  for(par in parameters){
    samples$conditional[[par]] <- .mix_samples(models, weights_list[[par]], converged, par, n_samples, seed)
  }
  for(par in c("mu", "sigma")){
    samples$averaged[[par]]    <- .mix_samples2(models, weights_all, converged, par, n_samples, seed)
  }
  for(par in c("mu", "sigma")[c("d", "r") %in% parameters]){
    samples$conditional[[par]] <- .mix_samples2(models, weights_list[[if(par == "mu") "d" else if(par == "sigma") "r"]], converged, par, n_samples, seed)
  }
  if(any(likelihoods %in% c("gamma0", "lognormal0"))){
    samples$averaged[["prop"]]    <- .mix_samples2(models, weights_all, converged, "prop", n_samples, seed)
    samples$conditional[["prop"]] <- .mix_samples2(models, weights_list[["p"]], converged, "prop", n_samples, seed)
  }


  prior_prob_all            <- prior_weights_all
  prior_prob_effect         <- sum(prior_weights_all[mm_d])
  prior_prob_heterogeneity  <- sum(prior_weights_all[mm_r])
  prior_prob_outliers       <- sum(prior_weights_all[mm_nu])
  prior_prob_proportion     <- sum(prior_weights_all[mm_p])
  
  posterior_prob_all            <- weights_all
  posterior_prob_effect         <- sum(weights_all[mm_d])
  posterior_prob_heterogeneity  <- sum(weights_all[mm_r])
  posterior_prob_outliers       <- sum(weights_all[mm_nu])
  posterior_prob_proportion     <- sum(weights_all[mm_p])
  
  
  # return the results
  output <- list(
    samples        = samples,
    BF             = list(
      effect           = BF_effect,
      heterogeneity    = BF_heterogeneity,
      outliers         = BF_outliers,
      proportion       = BF_proportion
    ),
    prior_prob     = list(
      all             = prior_prob_all,
      effect          = prior_prob_effect,
      heterogeneity   = prior_prob_heterogeneity,
      outliers        = prior_prob_outliers,
      proportion      = prior_prob_proportion
    ),
    posterior_prob = list(
      all             = posterior_prob_all,
      effect          = posterior_prob_effect,
      heterogeneity   = posterior_prob_heterogeneity,
      outliers        = posterior_prob_outliers,
      proportion      = posterior_prob_proportion
    )
  )
  return(output)
}
.mix_samples                <- function(models, weights, converged, parameter, n_samples, seed){
  
  if(is.null(weights)){
    return(NULL)
  }
  
  if(!is.null(seed)) set.seed(seed) else set.seed(1)
  samples <- NULL
  
  for(i in c(1:length(models))[converged]){

    # skip for missing nu parameter for all but t-distributions    
    if(round(n_samples * weights[i]) == 0 || is.null(models[[i]]$priors[[parameter]]))
      next
    
    if(models[[i]]$priors[[parameter]]$distribution == "point"){
      samples       <- c(samples, rep(models[[i]]$priors[[parameter]]$parameters$location, round(n_samples * weights[i])))
    }else{
      model_samples <- rstan::extract(models[[i]]$fit)
      ind           <- sample(nrow(model_samples[[parameter]]), round(n_samples * weights[i]), replace = TRUE)
      samples       <- c(samples, if(parameter == "nu") model_samples[[parameter]][ind] else model_samples[[parameter]][ind, 1])
    }
  }
  
  return(samples)
}
.mix_samples2               <- function(models, weights, converged, parameter, n_samples, seed){
  
  if(!is.null(seed)) set.seed(seed) else set.seed(1)
  samples <- NULL
  
  for(i in c(1:length(models))[converged]){
    
    if(round(n_samples * weights[i]) == 0)next
    
    model_samples <- rstan::extract(models[[i]]$fit)
    ind           <- sample(nrow(model_samples[[paste0(parameter, "_i")]]), round(n_samples * weights[i]), replace = TRUE)
    samples       <- rbind(samples, model_samples[[paste0(parameter, "_i")]][ind, ])
  }
  
  return(samples)
}
.compute_coeficients        <- function(RoBTT){
  return(c(
    "d"   = if(length(RoBTT$samples$averaged$d) != 0)  mean(RoBTT$samples$averaged$d),
    "r"   = if(length(RoBTT$samples$averaged$r) != 0)  mean(RoBTT$samples$averaged$r)
  ))
}
.inclusion_BF               <- function(prior_weights, posterior_weights, conditional_models){
  (sum(posterior_weights[conditional_models])/sum(posterior_weights[!conditional_models]))  /
    (sum(prior_weights[conditional_models])/sum(prior_weights[!conditional_models]))
}
.BF_format                  <- function(BF, BF01 = FALSE, logBF = FALSE){
  BF[is.nan(BF)] <- NA
  if(BF01){
    BF <- 1/BF
  }else{
    BF <- BF
  }
  if(logBF){
    BF <- log(BF)
  }
  return(BF)
}
.get_converged_models       <- function(object){

  # TODO
  return(rep(T, length(object$models)))
  
  converged <- NULL

  # basic convergence checks
  for(i in 1:length(object$models)){
    if(!.is_model_constant(object$models[[i]]$priors)){

      if(any(class(object$models[[i]]$fit) %in% c("simpleError", "error")) | is.infinite(object$models[[i]]$marg_lik$logml) | is.na(object$models[[i]]$marg_lik$logml)){
        converged <- c(converged, FALSE)
      }else{
        converged <- c(converged, TRUE)
      }

    }else{
      converged <- c(converged, TRUE)
    }
  }

  object$models <- object$models[converged]

  # remove models with unsatisfactory performance
  if(!is.null(object$control$allow_max_error) |!is.null(object$control$allow_max_rhat) | !is.null(object$control$allow_min_ESS)){
    diagnostics_summary <- summary.RoBTT(object, type = "models", diagnostics = TRUE, include_theta = object$control$allow_inc_theta)$diagnostics

    # deal with NAs for null models
    diagnostics_summary$"max(MCMC error)"[is.na(diagnostics_summary$"max(MCMC error)")] <- 0
    diagnostics_summary$"max(Rhat)"[is.na(diagnostics_summary$"max(Rhat)")]             <- 0
    diagnostics_summary$"min(ESS)"[is.na(diagnostics_summary$"min(ESS)")]               <- Inf


    if(!is.null(object$control$allow_max_error)){
      converged <- converged & (diagnostics_summary$"max(MCMC error)" < object$control$allow_max_error)
    }
    if(!is.null(object$control$allow_max_Rhat)){
      converged <- converged & diagnostics_summary$"max(Rhat)" < object$control$allow_max_rhat
    }
    if(!is.null(object$control$allow_min_ESS)){
      converged <- converged & diagnostics_summary$"min(ESS)"  > object$control$allow_min_ESS
    }
  }

  return(converged)
}
.balance_prob               <- function(object, converged_models){

  # extract data
  mm_mu      <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "mu"))
  mm_tau     <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "tau"))
  mm_omega   <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "omega"))
  prior_odds <- sapply(object$models, function(m)m$prior_odds_set)

  # check whether there is a comparable model for each non-converged models
  for(i in c(1:length(object$models))[!converged_models]){

    temp_ind  <- c(1:length(object$models))[-i]
    temp_same <- temp_ind[mm_mu[-i] == mm_mu[i] & mm_tau[-i] == mm_tau[i] & mm_omega[-i] == mm_omega[i] & converged_models[-i]]

    # if yes, transfer the prior odds
    if(length(temp_same) >= 1){
      prior_odds[temp_same] <- prior_odds[temp_same] + prior_odds[i] / length(temp_same)
      prior_odds[i] <- 0
      object$add_info$warnings <- c(object$add_info$warnings, "Some of the models failed to converge. However, there were other models with the same combination of presence/absence of effect/heterogeneity/publication bias and their prior probability was increased to account for the failed models.")
    }else{
      prior_odds[i] <- 0
      object$add_info$warnings <- c(object$add_info$warnings, "Some of the models failed to converge and their prior probability couldn't be balanced over models with the same combination of presence/absence of effect/heterogeneity/publication bias since they don't exist.")
    }
  }

  for(i in 1:length(object$models)){
    object$models[[i]]$prior_odds <- prior_odds[i]
  }

  return(object)
}
.model_refit_warnings       <- function(metadata){

  new_warn <- NULL

  # extract meta-data with fit-refit information
  refit_info <- t(sapply(metadata, function(x){
    if(is.null(x$refit_info)){
      return(c(x$i, NA))
    }else{
      return(c(x$i, x$refit_info))
    }
  }))

  marglik_info <- t(sapply(metadata, function(x){
    if(is.null(x$marg_lik)){
      return(c(x$i, NA))
    }else{
      return(c(x$i, x$marg_lik))
    }
  }))

  if(is.null(dim(refit_info)))  refit_info   <- matrix(refit_info,   ncol = 2)
  if(is.null(dim(marglik_info)))marglik_info <- matrix(marglik_info, ncol = 2)


  if(length(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Initial fit of %1$s %2$s failed due to incompatible starting values (most likely due to an outlier in the data and limited precision of t-distribution). Starting values for the mean parameter were therefore set to the mean of supplied data.",
      ifelse(length(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1]) == 1, "model", "models"),
      paste(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Initial fit of %1$s %2$s failed due to incompatible starting values (most likely due to an outlier in the data and limited precision of t-distribution). Starting values for the mean parameter were therefore set to the mean of supplied data and the model was refitted using boost likelihood function.",
      ifelse(length(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1]) == 1, "model", "models"),
      paste(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(refit_info[!refit_info[, 2] %in% c("empirical init", "refit with boost") & !is.na(refit_info[,2]), 1]) > 0){
    refit_info_messages_i <- refit_info[refit_info[, 2] != "not enough iterations" & !is.na(refit_info[,2]), 1]
    refit_info_messages   <- refit_info[refit_info[, 2] != "not enough iterations" & !is.na(refit_info[,2]), 2]

    for(i in 1:length(refit_info_messages_i)){
      new_warn <- c(new_warn, paste0("Model ", refit_info_messages_i[i]," failed with the following error: ", refit_info_messages[i]))
    }
  }


  if(length(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Marginal likelihood computation of %1$s %2$s couldn't be completed within the specified number of iterations.",
      ifelse(length(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1]) == 1, "model", "models"),
      paste(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 1]) > 0){
    marglik_info_messages_i <- marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 1]
    marglik_info_messages   <- marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 2]

    for(i in 1:length(marglik_info_messages_i)){
      new_warn <- c(new_warn, paste0("Marginal likelihood computation of model ", marglik_info_messages_i[i]," failed with the following error: ", marglik_info_messages[i]))
    }
  }

  return(new_warn)
}
.model_convergence_warnings <- function(object){

  new_warn <- NULL

  # used set values if specified by the user
  threshold_error <- ifelse(is.null(object$control$allow_max_error), Inf, object$control$allow_max_error)
  threshold_rhat  <- ifelse(is.null(object$control$allow_max_rhat), 1.05, object$control$allow_max_rhat)
  threshold_ESS   <- ifelse(is.null(object$control$allow_max_error), 100, object$control$allow_min_ESS)

  # get the diagnostics summary
  diagnostics_summary <- summary.RoBTT(object, type = "models", diagnostics = TRUE, include_theta = object$control$allow_inc_theta)$diagnostics

  # deal with NAs for null models
  diagnostics_summary$"max(MCMC error)"[is.na(diagnostics_summary$"max(MCMC error)")] <- 0
  diagnostics_summary$"max(Rhat)"[is.na(diagnostics_summary$"max(Rhat)")]             <- 0
  diagnostics_summary$"min(ESS)"[is.na(diagnostics_summary$"min(ESS)")]               <- Inf

  # find the problematic models
  warning_error <- rownames(diagnostics_summary)[diagnostics_summary$"max(MCMC error)" > threshold_error]
  warning_rhat  <- rownames(diagnostics_summary)[diagnostics_summary$"max(Rhat)"       > threshold_rhat]
  warning_ESS   <- rownames(diagnostics_summary)[diagnostics_summary$"min(ESS)"        < threshold_ESS]

  # add warnings messages
  if(length(warning_error) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with MCMC error larger than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_error) == 1, "Model", "Models"),
      paste(warning_error, collapse = ", "),
      threshold_error
    ))
  }

  if(length(warning_rhat) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with R-hat larger than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_rhat) == 1, "Model", "Models"),
      paste(warning_rhat, collapse = ", "),
      threshold_rhat
    ))
  }

  if(length(warning_ESS) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with ESS lower than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_ESS) == 1, "Model", "Models"),
      paste(warning_ESS, collapse = ", "),
      threshold_ESS
    ))
  }

  return(new_warn)
}

### helper functions for settings
.set_priors             <- function(prior_d, prior_r, prior_nu, prior_p, prior_d_null, prior_r_null, prior_nu_null, prior_p_null){

  priors     <- list()
  priors$d   <- .set_parameter_priors(prior_d_null,   prior_d,   "d")
  priors$r   <- .set_parameter_priors(prior_r_null,   prior_r,   "r")
  priors$nu  <- .set_parameter_priors(prior_nu_null,  prior_nu,  "nu")
  priors$p <- .set_parameter_priors(prior_p_null, prior_p, "p")

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
    if(class(priors_null) == "RoBTT.prior"){
      priors_null$is_null <- TRUE
    }else{
      stop(paste0("The null prior distribution for the ", parameter, " was not specified correctly."))
    }
  }
  if(is.null(priors_alt)){
    priors_alt <- NULL
  }else{
    # check that the prior is a prior object
    if(class(priors_alt) == "RoBTT.prior"){
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

  if(parameter %in% c("d", "p")){

    # check that the passed priors are supported for the parameter
    if(length(priors) > 0){
      for(i in 1:length(priors)){
        if(!priors[[i]]$distribution %in% c("normal", "lognormal", "t", "gamma", "invgamma", "point", "uniform", "beta", "exp"))
          stop(paste0(priors[[i]]$distribution," prior distribution is not supported for the ", parameter," parameter. See '?prior' for further information."))
      }
    }


  }else if(parameter == "r"){

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
  for(d in priors$d){
    for(r in priors$r){
      for(likelihood in likelihoods){
        if(likelihood == "t"){
          for(nu in priors$nu){
            models <- c(models, list(.create_model(d, r, nu, NULL, likelihood, d$prior_odds * r$prior_odds * nu$prior_odds)))
          }
        }else if(likelihood %in% c("gamma0", "lognormal0")){
          for(p in priors$p){
            models <- c(models, list(.create_model(d, r, NULL, p, likelihood, d$prior_odds * r$prior_odds * p$prior_odds)))
          }
        }else{
          models <- c(models, list(.create_model(d, r, NULL, NULL, likelihood, d$prior_odds * r$prior_odds)))
        } 
      }
    }
  }
  

  return(models)
}
.create_model           <- function(prior_d, prior_r, prior_nu, prior_p, likelihood, prior_odds){

  priors <- list()
  
  priors$d   <- prior_d
  priors$r   <- prior_r
  priors$nu  <- prior_nu
  priors$p <- prior_p

  # possibly simplify t to normal
  if(likelihood == "t" && prior_nu$distribution == "point" && prior_nu$parameters$location == Inf){
    likelihood <- "normal"
  }
  
  model <- list(
    priors         = priors,
    prior_odds     = prior_odds,
    prior_odds_set = prior_odds,
    likelihood     = likelihood
  )
  class(model) <- "RoBTT.model"

  return(model)
}
.set_control            <- function(control, chains, iter, warmup, thin, seed, parallel){

  # set the control list
  if(is.null(control)){
    control$bridge_max_iter <- 10000

    control$allow_max_error <- NULL
    control$allow_max_rhat  <- NULL
    control$allow_min_ESS   <- NULL

    control$silent          <- FALSE

    if(parallel){
      control$cores         <- parallel::detectCores() - 1
    }else{
      control$cores         <- 1
    }

  }else{
    if(is.null(control[["max_rhat"]])){
      control$max_rhat        <- 1.01
    }
    if(is.null(control[["bridge_max_iter"]])){
      control$bridge_max_iter <- 10000
    }
    if(is.null(control[["allow_max_error"]])){
      control$allow_max_error <- NULL
    }
    if(is.null(control[["allow_max_rhat"]])){
      control$allow_max_rhat  <- NULL
    }
    if(is.null(control[["allow_min_ESS"]])){
      control$allow_min_ESS   <- NULL
    }
    if(is.null(control[["silent"]])){
      control$silent          <- FALSE
    }
    if(is.null(control[["cores"]])){
      if(parallel){
        control$cores         <- parallel::detectCores() - 1
      }else{
        control$cores         <- 1
      }
    }

  }

  if(control[["cores"]] > 1){
    parallel <- TRUE
  }

  # add the main MCMC settings
  control$chains    <- chains
  control$iter      <- iter
  control$warmup    <- warmup
  control$thin      <- thin
  control$seed      <- seed
  control$parallel  <- parallel

  .check_control(control)
  return(control)
}
.check_control          <- function(control){
  # check whether only known controls were supplied
  known_controls <- c("chains", "iter", "warmup" , "adapt", "thin" ,"autofit", "max_error", "max_rhat", "max_time", "bridge_max_iter", "allow_max_error", "allow_max_rhat", "allow_min_ESS", "allow_inc_theta", "balance_prob", "silent", "progress_start", "progress_tick", "boost", "cores", "seed", "parallel", "effect_direction", "likelihood")
  if(any(!names(control) %in% known_controls))stop(paste0("The following control settings were not recognize: ", paste(names(control[!names(control) %in% known_controls]), collapse = ", ")))

  # check whether essential controls were supplied
  if(is.null(control[["chains"]])) stop("Number of chains must be defined.")
  if(is.null(control[["iter"]]))   stop("Number of iterations must be set.")
  if(is.null(control[["warmup"]])) stop("Number of warmup samples must be set.")
  if(is.null(control[["thin"]]))   stop("Thinning of the posterior samples must be set.")

  if(!is.numeric(control[["chains"]]) | !control[["chains"]] >= 1)  stop("At least one chains must be set.")
  if(!is.numeric(control[["iter"]])   | !control[["iter"]] >= 1)    stop("Number of iterations must be a positive number.")
  if(!is.numeric(control[["warmup"]]) | !control[["warmup"]] >= 1)  stop("Number of warmup samples must be a positive number.")
  if(!is.numeric(control[["thin"]])   | !control[["thin"]] >= 1)    stop("Thinning of the posterior samples must be a positive number.")
  if(!is.logical(control[["parallel"]]))                            stop("The usage of parallel evaluation must be a logical argument.")
  if(!is.numeric(control[["cores"]])  | !control[["cores"]] >= 1)   stop("Number of cores must be a positive number.")
  if(!is.numeric(control[["seed"]])   & !is.null(control[["seed"]]))stop("Seed must be a numeric argument.")


  # check convergence criteria
  if(!is.null(control[["allow_max_rhat"]]))  if(control[["allow_max_rhat"]] <= 1) stop("The maximum allowed R-hat must be higher than 1.")
  if(!is.null(control[["allow_min_ESS"]]))   if(control[["allow_min_ESS"]] <= 0)  stop("The minimum allowed ESS must be higher than 0.")

  if(control[["parallel"]]){
    if(!try(requireNamespace("parallel")))stop("parallel package needs to be installed for parallel processing. Run 'install.packages('parallel')'")
  }
  # now taken care of by the evaluation outside of R
  # runjags::runjags.options(silent.jags = control$silent, silent.runjags = control$silent)
}

# general helper functions
.is_parameter_null <- function(priors, par){
  return(if(is.null(priors[[par]])) TRUE else priors[[par]]$is_null)
}
.is_model_constant <- function(priors){

  constant <- NULL
  for(par in c("mu", "tau", "omega", "sigma", "p")){
    if(!is.null(priors[[par]])){
      constant <- c(constant, priors[[par]]$distribution == "point")
    }
  }

  constant <- all(constant)

  return(constant)
}
.get_no_support    <- function(models, par){

  no_support  <- NULL

  all_support <- sapply(models, function(m)m$priors[[par]]$truncation, simplify = FALSE)
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
