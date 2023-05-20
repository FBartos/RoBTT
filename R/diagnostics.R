#' @title Checks a fitted RoBTT object
#'
#' @description \code{diagnostics} creates visual
#' checks of individual models convergence. Numerical
#' overview of individual models can be obtained by
#' \code{summary(object, type = "models", diagnostics = TRUE)},
#' or even more detailed information by
#' \code{summary(object, type = "individual")}.
#'
#' @param fit a fitted RoBTT object
#' @param parameter a parameter to be plotted. Either
#' \code{"delta"}, \code{"rho"}, \code{"nu"}, \code{"mu"},
#' or \code{"sigma"}.
#' @param type type of MCMC diagnostic to be plotted.
#' Options are \code{"chains"} for the chains' trace plots,
#' \code{"autocorrelation"} for autocorrelation of the
#' chains, and \code{"densities"} for the overlaying
#' densities of the individual chains. Can be abbreviated to
#' first letters.
#' @param show_models MCMC diagnostics of which models should be
#' plotted. Defaults to \code{NULL} which plots MCMC diagnostics
#' for a specified parameter for every model that is part of the
#' ensemble.
#' @param title whether the model number should be displayed in title.
#' Defaults to \code{TRUE} when more than one model is selected.
#' @param lags number of lags to be shown for
#' \code{type = "autocorrelation"}. Defaults to \code{30}.
#' @param ... additional arguments to be passed to
#' \link[graphics]{par} if \code{plot_type = "base"}.
#'
#' @details The visualization functions are based on
#' \link[rstan]{stan_plot} function and its color schemes.
#'
#' @examples \dontrun{
#' # using the example data from Darwin
#' data("fertilization", package = "RoBTT")
#' fit <- RoBTT(
#'   x1       = fertilization$Self,
#'   x2       = fertilization$Crossed,
#'   prior_delta = prior("cauchy", list(0, 1/sqrt(2))),
#'   prior_rho   = prior("beta",   list(3, 3)),
#'   likelihood  = "normal",
#'   seed        = 1, 
#'   chains      = 1,
#'   warmup      = 1000,
#'   iter        = 2000,
#'   control     = set_control(adapt_delta = 0.95)
#' )
#'
#' ### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
#' # diagnostics function allows to visualize diagnostics of a fitted RoBTT object, for example,
#' # the trace plot for the mean parameter in each model model
#' diagnostics(fit, parameter = "delta", type = "chain")
#'
#' # in order to show the trace plot only for the 11th model, add show_models parameter
#' diagnostics(fit, parameter = "delta", type = "chain", show_models = 11)
#'
#' # furthermore, the autocorrelations
#' diagnostics(fit, parameter = "delta", type = "autocorrelation")
#'
#' # and overlying densities for each plot can also be visualize
#' diagnostics(fit, parameter = "delta", type = "densities")
#' }
#'
#'
#' @return \code{diagnostics} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object/list of objects (depending on the number of parameters to be plotted)
#' of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBTT()], [summary.RoBTT()]
#' @name diagnostics
#' @aliases diagnostics_autocorrelation diagnostics_trace diagnostics_density
#' @export diagnostics
#' @export diagnostics_density
#' @export diagnostics_autocorrelation
#' @export diagnostics_trace

#' @rdname diagnostics
diagnostics <- function(fit, parameter, type, show_models = NULL,
                        lags = 30, title = is.null(show_models) | length(show_models) > 1, ...){
  
  # check settings
  if(!is.RoBTT(fit))
    stop("Diagnostics are available only for RoBTT models.")
  BayesTools::check_char(parameter, "parameter")
  BayesTools::check_char(type, "type")
  
  # deal with bad type names
  if(substr(type, 1, 1) == "c"){
    type <- "trace"
  }else if(substr(type, 1, 1) == "t"){
    type <- "trace" # for trace
  }else if(substr(type, 1, 1) == "d"){
    type <- "density"
  }else if(substr(type, 1, 1) == "a"){
    type <- "autocorrelation"
  }else{
    stop("Unsupported diagnostics type. See '?diagnostics' for more details.")
  }
  
  if(parameter %in% c("delta", "rho", "nu", "mu")){
    parameter         <- parameter
    parameter_samples <- parameter
  }else if(parameter == "sigma"){
    parameter         <- "pooled_sigma"
    parameter_samples <- "pooled_sigma"
    
  }else{
    stop(paste0("The passed parameter does not correspond to any of main model parameter ('delta', 'rho', 'nu', 'mu', 'sigma'). See '?plot.RoBTT' for more details."))
  }
  
  # do the plotting
  models_ind <- 1:length(fit$models)
  if(!is.null(show_models)){
    models_ind <- models_ind[show_models]
  }
  
  plots <- list()
  
  for(i in models_ind){
    
    prior_names <- names(attr(fit$models[[i]][["fit"]], "prior_list"))
    prior_names <- prior_names[!sapply(attr(fit$models[[i]][["fit"]], "prior_list"), BayesTools::is.prior.point)]
    model_parameters <- c(prior_names, "mu", "pooled_sigma")
    
    if(!parameter_samples %in% model_parameters){
      
      plots[[i]] <- NULL
      
    }else if(inherits(fit$models[[i]][["fit"]], "null_model")){
      
      plots[[i]] <- NULL
      
    }else{

      if(type == "density"){
        plots[[i]] <- rstan::stan_dens(fit$models[[i]]$fit, pars = parameter) + 
          ggplot2::ylab("Density") + 
          ggplot2::xlab(.plot.RoBTT_par_names(parameter))
      }else if(type == "trace"){
        plots[[i]] <- rstan::traceplot(fit$models[[i]]$fit, pars = parameter) + 
          ggplot2::ylab(.plot.RoBTT_par_names(parameter)) + 
          ggplot2::xlab("Iteration")
      }else if(type == "autocorrelation"){
        plots[[i]] <- rstan::stan_ac(fit$models[[i]]$fit, pars = parameter, lags = lags, separate_chains = TRUE) + 
          ggplot2::ylab(substitute("Autocorrelation"~(x), list(x = .plot.RoBTT_par_names(parameter, quote = TRUE)))) +
          ggplot2::xlab("Lag")
      }
    }
  }
  
  
  # return the plots
  if(length(plots) == 0){
    plots <- NULL
  }else if(length(models_ind) == 1){
    plots <- plots[[models_ind]]
  }
  return(plots)
}


#' @rdname diagnostics
diagnostics_autocorrelation <- function(fit, parameter = NULL, plot_type = "base", show_models = NULL,
                                        lags = 30, title = is.null(show_models) | length(show_models) > 1, ...){
  diagnostics(fit = fit, parameter = parameter, type = "autocorrelation", plot_type = plot_type, show_models = show_models, lags = lags, title = title, ...)
}

#' @rdname diagnostics
diagnostics_trace           <- function(fit, parameter = NULL, plot_type = "base", show_models = NULL,
                                        title = is.null(show_models) | length(show_models) > 1, ...){
  diagnostics(fit = fit, parameter = parameter, type = "trace", plot_type = plot_type, show_models = show_models, title = title, ...)
}

#' @rdname diagnostics
diagnostics_density         <- function(fit, parameter = NULL, plot_type = "base", show_models = NULL,
                                        title = is.null(show_models) | length(show_models) > 1, ...){
  diagnostics(fit = fit, parameter = parameter, type = "density", plot_type = plot_type, show_models = show_models, title = title, ...)
}

