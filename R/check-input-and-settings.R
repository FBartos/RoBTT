#' @title Prints summary of \code{"RoBTT"} ensemble implied by the specified priors
#'
#' @description \code{check_setup} prints summary of \code{"RoBTT"} ensemble
#' implied by the specified prior distributions. It is useful for checking
#' the ensemble configuration prior to fitting all of the models.
#'
#' @inheritParams RoBTT
#' @param models should the models' details be printed.
#' @param silent do not print the results.
#'
#' @return \code{check_setup} invisibly returns list of summary tables.
#'
#' @seealso [RoBTT()]
#' @export
check_setup <- function(
    prior_delta  = prior(distribution = "cauchy",  parameters = list(location = 0, scale = sqrt(2)/2)),
    prior_rho    = prior(distribution = "beta",    parameters = list(alpha = 1, beta = 1)),
    prior_nu     = prior(distribution = "exp",     parameters = list(rate = 1)),
    
    prior_delta_null  = prior(distribution = "spike",  parameters = list(location = 0)),
    prior_rho_null    = prior(distribution = "spike",  parameters = list(location = 0.5)),
    prior_nu_null     = prior_none(),
    
    prior_mu = NULL, prior_sigma2 = NULL, truncation = NULL,
    
    models = FALSE, silent = FALSE){
  
  object <- list()
  object$priors      <- .set_priors(prior_delta, prior_rho, prior_nu, prior_delta_null, prior_rho_null, prior_nu_null, prior_mu, prior_sigma2, !is.null(truncation))
  object$models      <- .get_models(object$priors)
  
  ### model types overview
  effect        <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "delta"))
  heterogeneity <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "rho"))
  outliers      <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "nu"))
  
  # number of model types
  n_models    <- c(
    effect        = sum(effect),
    heterogeneity = sum(heterogeneity),
    outliers      = sum(outliers)
  )
  
  # extract model weights
  prior_weights   <- sapply(object$models, function(m) m$prior_weights)
  # standardize model weights
  prior_weights   <- prior_weights / sum(prior_weights)
  # conditional model weights
  models_prior <- c(
    effect         = sum(prior_weights[effect]),
    heterogeneity  = sum(prior_weights[heterogeneity]),
    outliers       = sum(prior_weights[outliers])
  )
  
  # create overview table
  components <- data.frame(
    "models"     = n_models,
    "prior_prob" = models_prior
  )
  rownames(components) <- c("Effect", "Heterogeneity", "Outliers")
  
  class(components)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(components))
  attr(components, "type")      <- c("n_models", "prior_prob")
  attr(components, "rownames")  <- TRUE
  attr(components, "n_models")  <- length(object$models)
  attr(components, "title")     <- "Components summary:"
  attr(components, "footnotes") <- NULL
  attr(components, "warnings")  <- NULL
  
  object$components <- components
  
  ### model details
  if(models){
    priors_effect        <- sapply(1:length(object$models), function(i) print(object$models[[i]]$priors$delta, silent = TRUE))
    priors_heterogeneity <- sapply(1:length(object$models), function(i) print(object$models[[i]]$priors$rho,   silent = TRUE))
    priors_outliers      <- sapply(1:length(object$models), function(i) {
      if(is.null(object$models[[i]]$priors$nu)){
        return("")
      }else{
        print(object$models[[i]]$priors$nu,    silent = TRUE)
      }
    })
    
    prior_weights    <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_weights)
    prior_prob       <- prior_weights / sum(prior_weights)
    model_likelihood <- sapply(object[["models"]], function(m) m[["likelihood"]])
    
    summary <- cbind.data.frame(
      "Model"         = 1:length(object$models),
      "Distribution"  = model_likelihood,
      "delta"         = priors_effect,
      "rho"           = priors_heterogeneity,
      "nu"            = priors_outliers,
      "prior_prob"    = prior_prob
    )
    class(summary)             <- c("BayesTools_table", "BayesTools_ensemble_inference", class(summary))
    attr(summary, "type")      <- c("integer", "string",  rep("prior", 3), "prior_prob")
    attr(summary, "rownames")  <- FALSE
    attr(summary, "title")     <- "Models overview:"
    attr(summary, "footnotes") <- NULL
    attr(summary, "warnings")  <- NULL
    
    object$summary <- summary
  }
  
  
  if(!silent){
    cat("Robust Bayesian t-test (set-up)\n")
    print(components, quote = FALSE, right = TRUE)
    
    if(models){
      cat("\n")
      print(summary, quote = FALSE, right = TRUE)
    }
  }
  
  return(invisible(object))
}


#' @title Convergence checks of the fitting process
#'
#' @description Set values for the convergence checks of the fitting process.
#'
#' @param max_Rhat maximum value of the R-hat diagnostic.
#' Defaults to \code{1.05}.
#' @param min_ESS minimum estimated sample size.
#' Defaults to \code{500}.
#' @param adapt_delta tuning parameter of HMC.
#' Defaults to \code{0.80}.
#' @param max_treedepth tuning parameter of HMC.
#' Defaults to \code{15}.
#' @param bridge_max_iter maximum number of iterations for the
#' \link[bridgesampling]{bridge_sampler} function. Defaults to \code{10000}
#'
#' 
#' @return \code{set_control} returns a list of control settings
#' and \code{set_convergence_checks} returns a list of convergence checks settings.
#'
#' @export set_control
#' @export set_convergence_checks
#' @name RoBTT_control
#' @aliases set_control, set_convergence_checks
NULL

#' @rdname RoBTT_control
set_convergence_checks  <- function(max_Rhat = 1.05, min_ESS = 500){
  
  BayesTools::check_real(max_Rhat, "max_Rhat", lower = 1)
  BayesTools::check_real(min_ESS, "min_ESS", lower = 0)
  
  convergence_checks <- list(
    max_Rhat            = max_Rhat,
    min_ESS             = min_ESS
  )
  
  return(convergence_checks)
}
#' @rdname RoBTT_control
set_control             <- function(adapt_delta = 0.80, max_treedepth = 15, bridge_max_iter = 1000){
  
  BayesTools::check_real(adapt_delta, "adapt_delta", lower = 0, upper = 1)
  BayesTools::check_int(max_treedepth, "max_treedepth", lower = 1)
  BayesTools::check_int(bridge_max_iter, "bridge_max_iter", lower = 1)
  
  control <- list(
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    bridge_max_iter = bridge_max_iter
  )
  
  return(control)
  
}


.check_and_list_convergence_checks <- function(convergence_checks){
  
  BayesTools::check_list(convergence_checks, "convergence_checks", check_names = c("max_Rhat", "min_ESS"))
  BayesTools::check_real(convergence_checks[["max_Rhat"]], "max_Rhat", lower = 1)
  BayesTools::check_real(convergence_checks[["min_ESS"]], "min_ESS", lower = 0)
  
  return(convergence_checks)
}
.stan_check_and_list_fit_settings  <- function(chains, warmup, iter, thin, parallel, cores, silent, seed, control, check_mins = list(chains = 1, warmup = 50, iter = 50, thin = 1), call = ""){
  
  BayesTools::check_int(chains, "chains",  lower = check_mins[["chains"]],  call = call)
  BayesTools::check_int(warmup, "warmup",  lower = check_mins[["warmup"]],  call = call)
  BayesTools::check_int(iter,   "iter",    lower = min(check_mins[["iter"]], warmup + 1), call = call)
  BayesTools::check_int(thin,   "thin",    lower = check_mins[["thin"]],    call = call)
  BayesTools::check_list(control, "control", check_names = c("adapt_delta", "max_treedepth", "bridge_max_iter"))
  
  BayesTools::check_bool(parallel, "parallel",                call = call)
  BayesTools::check_int(cores,     "cores", lower = 1,        call = call)
  BayesTools::check_bool(silent,   "silent",                  call = call)
  BayesTools::check_int(seed,      "seed", allow_NULL = TRUE, call = call)
  
  if(!parallel){
    cores <- 1
  }else if(cores > RoBTT.get_option("max_cores")){
    cores <- RoBTT.get_option("max_cores")
  }
  
  if(is.null(control[["adapt_delta"]])){
    control[["adapt_delta"]] <- 0.80
  }else{
    BayesTools::check_real(control[["adapt_delta"]], "adapt_delta", lower = 0, upper = 1)
  }
  if(is.null(control[["max_treedepth"]])){
    control[["max_treedepth"]] <- 15
  }else{
    BayesTools::check_int(control[["max_treedepth"]], "max_treedepth", lower = 1)
  }
  if(is.null(control[["bridge_max_iter"]])){
    control[["bridge_max_iter"]] <- 1000
  }else{
    BayesTools::check_int(control[["bridge_max_iter"]], "bridge_max_iter", lower = 1)
  }
  
  return(invisible(list(
    chains   = chains,
    warmup   = warmup,
    iter     = iter,
    thin     = thin,
    parallel = parallel,
    cores    = cores,
    silent   = silent,
    adapt_delta     = control[["adapt_delta"]],
    max_treedepth   = control[["max_treedepth"]],
    bridge_max_iter = control[["bridge_max_iter"]],
    seed     = seed
  )))
}

.update_fit_control        <- function(old_fit_control, chains, warmup, iter, thin, parallel, cores, silent, seed, control){
  
  if(is.null(chains)){
    chains <- old_fit_control[["chains"]]
  }
  if(is.null(warmup)){
    warmup  <- old_fit_control[["warmup"]]
  }
  if(is.null(iter)){
    iter <- old_fit_control[["iter"]]
  }
  if(is.null(thin)){
    thin  <- old_fit_control[["thin"]]
  }
  if(is.null(parallel)){
    parallel <- old_fit_control[["parallel"]]
  }
  if(is.null(silent)){
    silent <- old_fit_control[["silent"]]
  }
  if(is.null(seed)){
    seed   <- old_fit_control[["seed"]]
  }
  if(is.null(control)){
    control <- list(
      adapt_delta     = old_fit_control[["adapt_delta"]],
      max_treedepth   = old_fit_control[["max_treedepth"]],
      bridge_max_iter = old_fit_control[["bridge_max_iter"]]
    )
  }
  
  new_fit_control <- .stan_check_and_list_fit_settings(chains = chains, warmup = warmup, iter = iter, thin = thin, parallel = parallel, cores = chains, silent = silent, seed = seed, control = control)
  
  return(new_fit_control)
}
.update_convergence_checks <- function(old_convergence_checks, convergence_checks){
  
  if(!is.null(convergence_checks[["max_Rhat"]])){
    max_Rhat <- convergence_checks[["max_Rhat"]]
  }else{
    max_Rhat <- old_convergence_checks[["max_Rhat"]]
  }
  if(!is.null(convergence_checks[["min_ESS"]])){
    min_ESS <- convergence_checks[["min_ESS"]]
  }else{
    min_ESS <- old_convergence_checks[["min_ESS"]]
  }
  
  new_convergence_checks <- set_convergence_checks(max_Rhat = max_Rhat, min_ESS = min_ESS)
  new_convergence_checks <- .check_and_list_convergence_checks(new_convergence_checks)
  
  return(new_convergence_checks)
}

.check_data <- function(x1, x2, mean1, mean2, sd1, sd2, N1, N2, truncation){
  
  if(!is.null(mean1) & !is.null(mean2) &  !is.null(sd1) &  !is.null(sd2) &  !is.null(N1) &  !is.null(N2)){
    
    BayesTools::check_real(mean1, "mean1")
    BayesTools::check_real(mean2, "mean2")
    BayesTools::check_real(sd1, "sd1", lower = 0, allow_bound = FALSE)
    BayesTools::check_real(sd2, "sd2", lower = 0, allow_bound = FALSE)
    BayesTools::check_int(N1, "N1", lower = 1)
    BayesTools::check_int(N2, "N2", lower = 1)
    
    data <- list(
      x1    = NA,
      x2    = NA,
      mean1 = mean1,
      mean2 = mean2,
      sd1   = sd1,
      sd2   = sd2,
      N1    = N1,
      N2    = N2
    )
    attr(data, "summary") <- TRUE

  }else if(!is.null(x1) & !is.null(x2)){
    
    BayesTools::check_real(x1, "x1", check_length = 0)
    BayesTools::check_real(x2, "x2", check_length = 0)
    
    if(length(x1) < 1) stop("'x1' must contain at least one observation.")
    if(length(x2) < 2) stop("'x2' must contain at least one observation.")
    
    data <- list(
      x1    = x1,
      x2    = x2,
      mean1 = NA,
      mean2 = NA,
      sd1   = NA,
      sd2   = NA,
      N1    = NA,
      N2    = NA
    )
    attr(data, "summary") <- FALSE
    
    mean1 <- mean(x1)
    mean2 <- mean(x2)
    sd1   <- stats::sd(x1)
    sd2   <- stats::sd(x2)
    N1    <- length(x1)
    N2    <- length(x2)
    
  }else{
    stop("Insufficient data provided.")
  }
  
  # add truncation
  # - a single named integer sigma truncates both groups equally according to sd distance from the mean (common/grouped typed)
  # - a vector of length two specifies a common truncation range
  # - two vectors of length two specify different truncation ranges for each group
  if(is.null(truncation)){
    
    data[["is_trunc"]] <- 0
    data[["trunc1"]]   <- numeric()
    data[["trunc2"]]   <- numeric()
    
  }else{
    
    BayesTools::check_list(truncation, "truncation", check_names = c("sigma", "sigma1", "sigma2", "x", "x1", "x2"), allow_other = FALSE, all_objects = FALSE, check_length = 0)
    
    if("sigma" %in% names(truncation)){
      
      if(attr(data, "summary"))
        stop("Truncation by sigma is not supported for summary data.")
      
      BayesTools::check_real(truncation[["sigma"]], "truncation::sigma", lower = 0, allow_bound = FALSE)
      
      data[["is_trunc"]] <- 1
      data[["trunc1"]]   <- mean(c(x1, x2)) - truncation[["sigma"]] * sd(c(x1, x2))
      data[["trunc2"]]   <- mean(c(x1, x2)) - truncation[["sigma"]] * sd(c(x1, x2))
      
    }else if(all(c("sigma1", "sigma2") %in% names(truncation))){
      
      if(attr(data, "summary"))
        stop("Truncation by sigma is not supported for summary data.")
      
      BayesTools::check_real(truncation[["sigma1"]], "truncation::sigma1", lower = 0, allow_bound = FALSE)
      BayesTools::check_real(truncation[["sigma2"]], "truncation::sigma2", lower = 0, allow_bound = FALSE)
      
      data[["is_trunc"]] <- 1
      data[["trunc1"]]   <- mean(x1) - truncation[["sigma1"]] * sd(x1)
      data[["trunc2"]]   <- mean(x2) - truncation[["sigma2"]] * sd(x2)
      
    }else if("x" %in% names(truncation)){
      
      BayesTools::check_real(truncation[["x"]], "truncation::x", check_length = 2)
      BayesTools::check_real(truncation[["x"]][1], "truncation::x[1]", upper = truncation[["x"]][2], allow_bound = FALSE)
      
      data[["is_trunc"]] <- 1
      data[["trunc1"]]   <- truncation[["x"]]
      data[["trunc2"]]   <- truncation[["x"]]
      
    }else if(all(c("x1", "x2") %in% names(truncation))){
      
      if(attr(data, "summary"))
        stop("Truncation by x is not supported for summary data.")
      
      BayesTools::check_real(truncation[["x1"]], "truncation::x1", check_length = 2)
      BayesTools::check_real(truncation[["x2"]], "truncation::x2", check_length = 2)
      BayesTools::check_real(truncation[["x1"]][1], "truncation::x1[1]", upper = truncation[["x1"]][2], allow_bound = FALSE)
      BayesTools::check_real(truncation[["x2"]][1], "truncation::x2[1]", upper = truncation[["x2"]][2], allow_bound = FALSE)
      
      data[["is_trunc"]] <- 1
      data[["trunc1"]]   <- truncation[["x1"]]
      data[["trunc2"]]   <- truncation[["x2"]]
      
    }
    
    # remove the truncated values
    if(!attr(data, "summary")){
      
      x1_outside <- x1 < data[["trunc1"]][1] | x1 > data[["trunc1"]][2]
      x2_outside <- x2 < data[["trunc2"]][1] | x2 > data[["trunc2"]][2]
      
      if(sum(x1_outside) > 0 | sum(x2_outside) > 0){
        warning(paste0("Truncation removed ", sum(x1_outside), " and ", sum(x2_outside), " observations from group 1 and 2, respectively."), immediate. = TRUE, call. = FALSE)
      }
      
      attr(data, "n_truncated") <- sum(x1_outside) + sum(x2_outside)
      
      if(sum(!x1_outside) < 1) stop("'x1' must contain at least one observation.")
      if(sum(!x2_outside) < 2) stop("'x2' must contain at least one observation.")
      
      data[["x1"]] <- x1[!x1_outside]
      data[["x2"]] <- x2[!x2_outside]
    }
    
  }
  
  return(data)
}

# some functions for the JASP implementation
.RoBTT_collect_dots      <- function(...){
  
  dots <- list(...)
  
  known_dots <- c("is_JASP")
  if(any(!names(dots) %in% known_dots))
    stop(paste0("Uknown arguments to 'RoBTT': ", paste("'", names(dots)[!names(dots) %in% known_dots], "'", collapse = ", "), "."), call. = FALSE)
  
  if(is.null(dots[["is_JASP"]])){
    dots[["is_JASP"]] <- FALSE
  }else{
    dots[["is_JASP"]] <- TRUE
  }
  
  return(dots)
}
.JASP_progress_bar_start <- function(n){
  eval(expr = parse(text = 'startProgressbar(n)'))
}
.JASP_progress_bar_tick  <- function(){
  eval(expr = parse(text = 'progressbarTick()'))
}
