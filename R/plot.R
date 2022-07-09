#' @title Plots a fitted RoBTT object
#'
#' @description \code{plot.RoBTT} allows to visualize
#' different \code{"RoBTT"} object parameters in various
#' ways. See \code{type} for the different model types.
#'
#' @param x a fitted RoBTT object
#' @param parameter a parameter to be plotted. Defaults to
#' \code{"delta"} (for the effect size). The additional options
#' are \code{"rho"} (for the heterogeneity),
#' \code{"nu"} (for the degrees of freedom).
#' @param conditional whether conditional estimates should be
#' plotted. Defaults to \code{FALSE} which plots the model-averaged
#' estimates. Note that both \code{"weightfunction"} and
#' \code{"PET-PEESE"} are always ignoring the other type of
#' publication bias adjustment.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, and y-axis range
#' respectively.
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' }
#'
#'
#' @return \code{plot.RoBTT} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBTT()]
#' @export
plot.RoBTT  <- function(x, parameter = "mu",
                        conditional = FALSE, plot_type = "base", prior = FALSE, dots_prior = NULL, ...){
  
  # check whether plotting is possible
  if(!all(.get_model_convergence(x)))
    stop("Some of the models did not converge.")
  
  #check settings
  BayesTools::check_char(parameter, "parameter")
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_char(parameter, "plot_type", allow_values = c("delta", "rho", "nu"))
  
  # choose the samples
  if(conditional){
    samples <- x[["RoBTT"]][["posteriors_conditional"]]
  }else{
    samples <- x[["RoBTT"]][["posteriors"]]
  }
  
  if(conditional && is.null(samples[[parameter]])){
    switch(
      parameter,
      "delta" = stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the effect. Please, verify that you specified at least one model assuming the presence of the effect."),
      "rho"   = stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the heterogeneity. Please, verify that you specified at least one model assuming the presence of the heterogeneity."),
      "nu"    = stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of outliers. Please, verify that you specified at least one model assuming the presence of outliers.")
    )
  }else if(is.null(samples[[parameter]])){
    switch(
      parameter,
      "delta" = stop("The ensemble does not contain any posterior samples model-averaged across the effect. Please, verify that you specified at least one model for the effect."),
      "rho"   = stop("The ensemble does not contain any posterior samples model-averaged across the heterogeneity. Please, verify that you specified at least one model for the heterogeneity."),
      "nu"    = stop("The ensemble does not contain any posterior samples model-averaged across outliers. Please, verify that you specified at least one model for outliers.")
    )
  }
  
  # pretend that infinite degrees of freedom are 0 to make plotting possible for nu
  if(parameter == "nu"){
    
    if(prior)
      stop("Prior and posterior plots are not implemented for the degrees of freedom parameter.")
    # the prior needs to be "shifted by two to the right"
    
    samples[["nu"]][is.infinite(samples[["nu"]])] <- 0
    
    for(i in seq_along(x[["models"]])){
      if(x[["models"]][[i]][["likelihood"]] == "normal"){
        attr(samples[["nu"]], "prior_list")[[i]] <- prior("spike", list("location" = 0))
      }
    }
    
  }
  
  dots       <- .set_dots_plot(...)
  dots_prior <- .set_dots_prior(dots_prior)
  
  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- parameter
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$transformation           <- NULL
  args$transformation_arguments <- NULL
  args$transformation_settings  <- FALSE
  args$rescale_x                <- FALSE
  args$par_name                 <- .plot.RoBTT_par_names(parameter, x)
  args$dots_prior               <- dots_prior
  
  plot <- do.call(BayesTools::plot_posterior, args)
  
  
  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}


.set_dots_plot        <- function(...){
  
  dots <- list(...)
  if(is.null(dots[["col"]])){
    dots[["col"]]      <- "black"
  }
  if(is.null(dots[["col.fill"]])){
    dots[["col.fill"]] <- "#4D4D4D4C" # scales::alpha("grey30", .30)
  }
  
  return(dots)
}
.set_dots_prior       <- function(dots_prior){
  
  if(is.null(dots_prior)){
    dots_prior <- list()
  }
  
  if(is.null(dots_prior[["col"]])){
    dots_prior[["col"]]      <- "grey60"
  }
  if(is.null(dots_prior[["lty"]])){
    dots_prior[["lty"]]      <- 1
  }
  if(is.null(dots_prior[["col.fill"]])){
    dots_prior[["col.fill"]] <- "#B3B3B34C" # scales::alpha("grey70", .30)
  }
  
  return(dots_prior)
}
.plot.RoBTT_par_names <- function(par, fit){
  
  return(switch(
    par,
    "delta" = expression(delta),
    "rho"   = expression(rho),
    "nu"    = expression(nu)
  ))
}
