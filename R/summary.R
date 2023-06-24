#' @title Prints a fitted 'RoBTT' object
#'
#' @param x a fitted 'RoBTT' object.
#' @param ... additional arguments.
#'
#'
#' @return \code{print.RoBTT} invisibly returns the print statement.
#'
#' @seealso [RoBTT()]
#' @export
print.RoBTT <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEstimates:\n")
  print(stats::coef(x))
}


#' @title Summarize fitted 'RoBTT' object
#'
#' @description \code{summary.RoBTT} creates summary tables for a
#' RoBTT object.
#'
#' @param object a fitted 'RoBTT' object
#' @param type whether to show the overall 'RoBTT' results (\code{"ensemble"}),
#' an overview of the individual models (\code{"models"}), an overview of
#' the individual models MCMC diagnostics (\code{"diagnostics"}), or a detailed summary
#' of the individual models (\code{"individual"}). Can be abbreviated to first letters.
#' @param conditional show the conditional estimates (assuming that the
#' alternative is true). Defaults to \code{FALSE}. Only available for
#' \code{type == "conditional"}.
#' @param group_estimates show the model-averaged mean and standard deviation estimates for each group.
#' @param probs quantiles of the posterior samples to be displayed.
#' Defaults to \code{c(.025, .975)}
#' @param logBF show log of Bayes factors. Defaults to \code{FALSE}.
#' @param BF01 show Bayes factors in support of the null hypotheses. Defaults to
#' \code{FALSE}.
#' @param short_name whether priors names should be shortened to the first
#' (couple) of letters. Defaults to \code{FALSE}.
#' @param remove_spike_0 whether spike prior distributions with location at zero should
#' be omitted from the summary. Defaults to \code{FALSE}.
#' @param ... additional arguments
#'
#' @examples
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
#'
#' # estimates from the conditional models can be obtained with
#' summary(fit, conditional = TRUE)
#'
#' # overview of the models and their prior and posterior probability, marginal likelihood,
#' # and inclusion Bayes factor can be obtained with
#' summary(fit, type = "models")
#'
#' # diagnostics overview, containing the maximum R-hat, minimum ESS, maximum MCMC error, and
#' # maximum MCMC error / sd across parameters for each individual model can be obtained with
#' summary(fit, type = "diagnostics")
#'
#' # summary of individual models and their parameters can be further obtained by
#' summary(fit, type = "individual")
#'
#'
#' @return \code{summary.RoBTT} returns a list of tables of class 'BayesTools_table'.
#'
#' @seealso [RoBTT()]
#' @export
summary.RoBTT       <- function(object, type = "ensemble", conditional = FALSE,
                                group_estimates = FALSE, probs = c(.025, .975), logBF = FALSE, BF01 = FALSE,
                                short_name = FALSE, remove_spike_0 = FALSE, ...){
  
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(type, "type")
  BayesTools::check_bool(group_estimates, "group_estimates")
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_bool(BF01,  "BF01")
  BayesTools::check_bool(logBF, "logBF")
  BayesTools::check_bool(short_name, "short_name")
  BayesTools::check_bool(remove_spike_0, "remove_spike_0")
  
  # print diagnostics if all models fail to converge
  if(any(!.get_model_convergence(object))){
    if(substr(type,1,1) != "d")
      warning("At least one model failed to converge. Model diagnostics were printed instead.")
    type        <- "diagnostics"
  }
  
  if(substr(type,1,1) == "e"){
    
    # obtain components overview
    components <- BayesTools::ensemble_inference_table(
      inference  = object$RoBTT[["inference"]],
      parameters = names(object$RoBTT[["inference"]]),
      logBF      = logBF,
      BF01       = BF01,
      title      = "Components summary:"
    )
    
    # obtain estimates tables
    estimates <- BayesTools::ensemble_estimates_table(
      samples    = object$RoBTT[["posteriors"]],
      parameters = names(object$RoBTT[["posteriors"]]),
      probs      = probs,
      title      = "Model-averaged estimates:",
      warnings   = .collect_errors_and_warnings(object)
    )
    
    if(group_estimates){
      estimates_group <- BayesTools::ensemble_estimates_table(
        samples    = object$RoBTT[["posteriors_est"]],
        parameters = names(object$RoBTT[["posteriors_est"]]),
        probs      = probs,
        title      = "Model-averaged group parameter estimates:",
        warnings   = .collect_errors_and_warnings(object)
      )
    }
    
    # deal with possibly empty table in case of no alternative models
    if(is.null(object$RoBTT[["posteriors_conditional"]])){
      estimates_conditional                    <- data.frame(matrix(nrow = 0, ncol = length(probs) + 2))
      colnames(estimates_conditional)          <- c("Mean", "Median", probs)
      class(estimates_conditional)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(estimates_conditional))
      attr(estimates_conditional, "type")      <- rep("estimate", ncol(estimates_conditional))
      attr(estimates_conditional, "rownames")  <- TRUE
      attr(estimates_conditional, "title")     <- "Conditional estimates:"
      attr(estimates_conditional, "warnings")  <- .collect_errors_and_warnings(object)
    }else{
      estimates_conditional <- BayesTools::ensemble_estimates_table(
        samples    = object$RoBTT[["posteriors_conditional"]],
        parameters = names(object$RoBTT[["posteriors_conditional"]]),
        probs      = probs,
        title      = "Conditional estimates:",
        warnings   = .collect_errors_and_warnings(object)
      )
    }
    
    
    ### return results
    output <- list(
      call       = object[["call"]],
      title      = "Robust Bayesian t-test",
      components = components,
      estimates  = estimates
    )
    
    if(group_estimates){
      output$estimates_group <- estimates_group
    }
    
    if(conditional){
      output$estimates_conditional <- estimates_conditional
    }
    
    class(output) <- "summary.RoBTT"
    attr(output, "type") <- "ensemble"
    
    return(output)
    
  }else if(substr(type,1,1) == "m"){
    
    components <- names(object$RoBTT[["inference"]])
    parameters <- names(object$RoBTT[["posteriors"]])
    
    summary <- BayesTools::ensemble_summary_table(
      models         = object[["models"]],
      parameters     = parameters,
      title          = "Models overview:",
      footnotes      = NULL,
      warnings       = .collect_errors_and_warnings(object),
      short_name     = short_name,
      remove_spike_0 = remove_spike_0
    )
    
    summary <- BayesTools::add_column(
      summary,
      column_title    = "Distribution",
      column_values   = sapply(object[["models"]], function(m) m[["likelihood"]]),
      column_position = 2,
      column_type     = "string")
    
    output <- list(
      call       = object[["call"]],
      title      = "Robust Bayesian t-test",
      summary    = summary
    )
    
    class(output) <- "summary.RoBTT"
    attr(output, "type") <- "models"
    
    return(output)
    
  }else if(substr(type,1,1) == "d"){
    
    components <- names(object$RoBTT[["inference"]])
    parameters <- names(object$RoBTT[["posteriors"]])
    
    diagnostics <- BayesTools::ensemble_diagnostics_table(
      models         = object[["models"]],
      parameters     = parameters,
      title          = "Diagnostics overview:",
      footnotes      = NULL,
      warnings       = .collect_errors_and_warnings(object),
      short_name     = short_name,
      remove_spike_0 = remove_spike_0
    )
    
    diagnostics <- BayesTools::add_column(
      diagnostics,
      column_title    = "Distribution",
      column_values   = sapply(object[["models"]], function(m) m[["likelihood"]]),
      column_position = 2,
      column_type     = "string")
    
    output <- list(
      call        = object[["call"]],
      title       = "Robust Bayesian t-test",
      diagnostics = diagnostics
    )
    
    class(output) <- "summary.RoBTT"
    attr(output, "type") <- "diagnostics"
    
    return(output)
    
  }else if(substr(type, 1, 1) == "i"){
    
    output <- list(
      call       = object[["call"]],
      title      = "Robust Bayesian t-test",
      models     = list()
    )
    
    for(i in seq_along(object[["models"]])){
      
      summary  <- BayesTools::model_summary_table(
        model          = object[["models"]][[i]],
        short_name     = short_name,
        remove_spike_0 = remove_spike_0
      )

      estimates <- object[["models"]][[i]][["fit_summary"]]
      attr(estimates, "warnings")  <- object[["models"]][[i]][["warnings"]]
      attr(estimates, "title")     <- "Parameter estimates:"
      
      output[["models"]][[i]] <- list(
        summary   = summary,
        estimates = estimates
      )
    }
    
    class(output) <- "summary.RoBTT"
    attr(output, "type") <- "individual"
    
    return(output)
    
  }else{
    stop(paste0("Unknown summary type: '", type, "'."))
  }
}


#' @title Prints summary object for 'RoBTT' method
#'
#' @param x a summary of a 'RoBTT' object
#' @param ... additional arguments
#'
#'
#' @return \code{print.summary.RoBTT} invisibly returns the print statement.
#'
#' @seealso [RoBTT()]
#' @export
print.summary.RoBTT <- function(x, ...){
  
  cat("Call:\n")
  print(x[["call"]])
  
  cat("\n")
  cat(x[["title"]])
  
  
  if(attr(x, "type") == "ensemble"){
    
    cat("\n")
    print(x[["components"]])
    
    cat("\n")
    print(x[["estimates"]])
    
    if(!is.null(x[["estimates_group"]])){
      cat("\n")
      print(x[["estimates_group"]])
    }
    
    if(!is.null(x[["estimates_conditional"]])){
      cat("\n")
      print(x[["estimates_conditional"]])
    }
    
    return(invisible())
    
  }else if(attr(x, "type") == "models"){
    
    cat("\n")
    print(x[["summary"]])
    
    return(invisible())
    
  }else if(attr(x, "type") == "diagnostics"){
    
    cat("\n")
    print(x[["diagnostics"]])
    
    return(invisible())
    
  }else if(attr(x, "type") == "individual"){
    
    for(i in seq_along(x[["models"]])){
      
      if(i > 1){
        cat("\n")
      }
      print(x[["models"]][[i]][["summary"]])
      
      cat("\n")
      print(x[["models"]][[i]][["estimates"]])
    }
    
    return(invisible())
  }
}


#' @title Reports whether x is a 'RoBTT' object
#'
#' @param x an object to test
#'
#'
#' @return \code{is.RoBTT} returns a boolean.
#'
#' @export
is.RoBTT            <- function(x){
  inherits(x, "RoBTT")
}



#' @title Interprets results of a 'RoBTT' model.
#'
#' @description \code{interpret} creates a brief textual summary
#' of a fitted 'RoBTT' object.
#'
#' @inheritParams summary.RoBTT
#'
#'
#' @return \code{interpret} returns a character.
#'
#' @export
interpret           <- function(object){
  
  specification <- list(
    list(
      inference           = "Effect",
      samples             = "delta",
      inference_name      = "effect",
      inference_BF_name   = "BF_10",
      samples_name        = "delta"
    ),
    list(
      inference           = "Heterogeneity",
      samples             = "rho",
      inference_name      = "heterogeneity",
      inference_BF_name   = "BF^rho",
      samples_name        = "rho"
    ),
    list(
      inference           = "Outliers",
      samples             = "nu",
      inference_name      = "outliers",
      inference_BF_name   = "BF^nu",
      samples_name        = "nu"
    )[]
  )
  specification <- specification[c(
    any(names(object$RoBTT[["inference"]]) == "Effect"),
    any(names(object$RoBTT[["inference"]]) == "Heterogeneity"),
    any(names(object$RoBTT[["inference"]]) == "Outliers")
  )]
  
  text <- BayesTools::interpret(
    inference     = object$RoBTT[["inference"]],
    samples       = object$RoBTT[["posteriors"]],
    specification = specification,
    method        = "Robust Bayesian t-test"
  )
  
  return(text)
}
