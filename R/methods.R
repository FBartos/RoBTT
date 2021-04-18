#' @title Prints a fitted RoBTT object
#'
#' @param x a fitted RoBTT object.
#' @param ... additional arguments.
#' @export  print.RoBTT
#' @rawNamespace S3method(print, RoBTT)
#' @seealso [RoBTT()]
print.RoBTT <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEstimates:\n")
  print(stats::coef(x))
}


#' @title Summarize fitted RoBTT object
#'
#' @description \code{summary.RoBTT} creates a numerical
#' summary of the RoBTT object.
#'
#' @param object a fitted RoBTT object.
#' @param type whether to show the overall RoBTT results (\code{"ensemble"}),
#' an overview of the individual models (\code{"models"}), or detailed summary
#' for the individual models (\code{"individual"}).
#' @param conditional show the conditional estimates (assuming that the
#' alternative is true). Defaults to \code{FALSE}. Only available for
#' \code{type == "conditional"}.
#' @param diagnostics show the maximum R-hat and minimum ESS for the main
#' parameters in each of the models. Only available for \code{type = "ensemble"}.
#' @param include_theta whether the estimated random effects should be included
#' either in the summaries.
#' @param probs quantiles of the posterior samples to be displayed.
#' Defaults to \code{c(.025, .50, .975)}
#' @param logBF show log of the BFs. Defaults to \code{FALSE}.
#' @param BF01 show BF in support of the null hypotheses. Defaults to
#' \code{FALSE}.
#' @param digits_estimates a number of decimals for rounding the estimates.
#' Defaults to \code{3}.
#' @param digits_BF a number of decimals for rounding the BFs. Defaults to \code{3}.
#' @param ... additional arguments
#'
#' @return summary of a RoBTT object
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBTT(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' # summary can provide many details about the model
#' summary(fit)
#'
#' # note that the summary function contains additional arguments
#' # that allow to obtain a specific output, i.e, the conditional estimates
#' # (assuming that the non-null models are true) can be obtained
#' summary(fit, conditional = TRUE)
#'
#' # overview of the models and their prior and posterior probability, marginal likelihood,
#' # and inclusion Bayes factor:
#' summary(fit, type = "models")
#'
#' # and the model diagnostics overview, containing maximum R-hat and minimum ESS across parameters
#' # but see '?diagnostics' for diagnostics plots for individual model parameters
#' summary(fit, type = "models", diagnostics = TRUE)
#'
#' # summary of individual models and their parameters can be further obtained by
#' summary(fit, type = "individual")
#'
#' }
#' @note See [diagnostics()] for visual convergence checks of the individual models.
#' @method summary RoBTT
#' @export summary.RoBTT
#' @rawNamespace S3method(summary, RoBTT)
#' @seealso [RoBTT()] [diagnostics()]
summary.RoBTT       <- function(object, type = if(diagnostics) "models" else "ensemble",
                                conditional = FALSE, diagnostics = FALSE,
                                probs = c(.025, .975), logBF = FALSE, BF01 = FALSE,
                                digits_estimates = 3, digits_BF = 3,...){


  # print diagnostics if all models fail to converge
  if(!any(object$add_info$converged)){
    if(substr(type,1,1) != "m" & !diagnostics)warning("All models failed to converge. Model diagnostics were printed instead.")
    type        <- "models"
    diagnostics <- TRUE
  }


  if(substr(type,1,1) == "e"){

    ### model estimates
    parameters <- list()
    estimates  <- list()

    mm_d  <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "d"))
    mm_r  <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "r"))
    mm_nu <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "nu"))
    
    pars <- c("d", "r", "nu")[c(sum(mm_d), sum(mm_r), sum(mm_nu))]
    
    # compute quantiles
    if(!is.null(probs)){
      if(length(probs) != 0){

        if(!is.numeric(probs) | !is.vector(probs))stop("The passed probabilities 'probs' must be a numeric vector.")
        if(!(all(probs > 0) & all(probs < 1)))stop("The passed probabilities 'probs' must be higher than 0 and lower than 1.")

        # quantiles
        for(type in c("averaged", "conditional")){
          
          parameters[[type]] <- data.frame(do.call(rbind, lapply(object$RoBTT$samples[[type]][pars], stats::quantile, probs = probs)))
          estimates[[type]]  <- data.frame(do.call(rbind, lapply(object$RoBTT$samples[[type]][c("mu", "sigma")],  function(x){
            t(apply(x, 2, stats::quantile, probs = probs))
          })))
          
          colnames(parameters[[type]]) <- probs
          colnames(estimates[[type]])  <- probs

        }
      }
    }

    # point estimates
    for(type in c("averaged", "conditional")){
      parameters[[type]] <- cbind(do.call(rbind, lapply(object$RoBTT$samples[[type]][pars], function(x)c("Mean" = mean(x), "Median" = stats::median(x)))), parameters[[type]])
      estimates[[type]]  <- cbind(do.call(rbind, lapply(object$RoBTT$samples[[type]][c("mu", "sigma")],  function(x){
        t(apply(x, 2, function(xi)c("Mean" = mean(xi), "Median" = stats::median(xi))))
      })), estimates[[type]])
    }


    ### fixing naming
    rownames(estimates$averaged)    <- c("mu[1]", "mu[2]", "sigma[1]", "sigma[2]")
    rownames(estimates$conditional) <- c("mu[1]", "mu[2]", "sigma[1]", "sigma[2]")


    ### model types overview
    parameters_models         <- sapply(pars, function(par){
      sum(sapply(object$models, function(m)!.is_parameter_null(m$priors, par)))
    })
    parameters_prior_prob     <- unlist(object$RoBTT$prior_prob[c("effect", "heterogeneity", "outliers")])
    parameters_posterior_prob <- unlist(object$RoBTT$posterior_prob[c("effect", "heterogeneity", "outliers")])
    parameters_BF             <- unlist(object$RoBTT$BF[c("effect", "heterogeneity", "outliers")])
    parameters_BF             <- .BF_format(parameters_BF, BF01, logBF)
    overview_tab              <- data.frame(parameters_models, parameters_prior_prob, parameters_posterior_prob, parameters_BF)
    colnames(overview_tab) <- c(
      "Models",
      "Prior prob.",
      "Post. prob.",
      paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")")
    )
    rownames(overview_tab) <- c("Effect", "Heterogeneity", "Outliers")
    
    if(!conditional){
      parameters[["conditional"]] <- NULL
      estimates[["conditional"]]  <- NULL
    }

    ### return results
    res <- list(
      call       = object$call,
      overview   = overview_tab,
      parameters = parameters,
      estimates  = estimates,
      add_info   = list(
        n_models         = length(object$models),
        digits_estimates = digits_estimates,
        digits_BF        = digits_BF,
        type             = "ensemble",
        failed           = sum(!object$add_info$converged)
      )
    )


  }else if(substr(type,1,1) == "m"){


    priors_d   <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors[["d"]],  silent = TRUE))
    priors_r   <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors[["r"]],  silent = TRUE))
    priors_nu  <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors[["nu"]], silent = TRUE))


    prior_odds     <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_odds)
    prior_prob     <- prior_odds / sum(prior_odds)
    marg_lik       <- sapply(1:length(object$models), function(i)object$models[[i]]$marg_lik$logml)
    posterior_prob <- bridgesampling::post_prob(marg_lik, prior_prob = prior_prob)
    BF             <- sapply(1:length(object$models), function(i){
      temp_mm <- rep(FALSE, length(object$models))
      temp_mm[i] <- TRUE
      .inclusion_BF(prior_prob, posterior_prob, temp_mm)
    })

    BF[is.nan(BF)] <- NA
    if(BF01){
      BF <- 1/BF
    }else{
      BF <- BF
    }
    if(logBF){
      BF <- log(BF)
    }

    overview_tab <- data.frame(
      priors_d,
      priors_r,
      priors_nu,
      prior_prob,
      posterior_prob,
      marg_lik,
      BF,
      stringsAsFactors = FALSE
    )
    rownames(overview_tab) <- NULL
    colnames(overview_tab) <- c("Prior d", "Prior r", "Prior nu", "Prior prob.", "Post. prob.", "log(MargLik)",
                                paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")"))

    # add the summary model diagnostics
    if(diagnostics){

      diagnostics_tab <- overview_tab[,1:3]

      # extract max(R-hat) & min(ESS)
      diag_sum <- sapply(1:length(object$models), function(i){

        temp_x <- object$models[[i]]$fit

        if(length(temp_x) == 0 | any(class(object$models[[i]]$fit) %in% c("simpleError","error"))){
          return(c(NA, NA, NA))
        }else{
          s.x <- object$models[[i]]$fit_summary
          return(c(
            Rhat  = max(s.x[,"Rhat"]),
            ESS   = min(s.x[,"n_eff"])
          ))
        }
       })

      diagnostics_tab$"max(Rhat)"        <- diag_sum[1,]
      diagnostics_tab$"min(ESS)"         <- diag_sum[2,]

    }


    res <- list(
      call     = object$call,
      overview = overview_tab,
      add_info = list(
        n_models         = length(object$models),
        digits_estimates = digits_estimates,
        digits_BF        = digits_BF,
        type             = "models"
      )
    )

    if(diagnostics){
      res$diagnostics <- diagnostics_tab
    }


  }else if(substr(type, 1, 1) == "i"){

    overview_tabs <- list()

    prior_odds     <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_odds)
    prior_prob     <- prior_odds / sum(prior_odds)
    posterior_prob <- object$RoBTT$posterior_prob$all
    marg_lik       <- sapply(1:length(object$models), function(i)object$models[[i]]$marg_lik$logml)
    BF             <- sapply(1:length(object$models), function(i){
      temp_mm <- rep(FALSE, length(object$models))
      temp_mm[i] <- TRUE
      .inclusion_BF(prior_prob, posterior_prob, temp_mm)
    })
    BF             <- .BF_format(BF, BF01, logBF)


    for(i in 1:length(object$models)){

      if(length(object$models[[i]]$fit) == 0 |  any(class(object$models[[i]]$fit) %in% c("simpleError","error"))){

        s.x <- NULL

      }else{

        s.x <- object$models[[i]]$fit_summary[,c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat"), drop = FALSE]

        colnames(s.x) <- c("Mean", "SD", ".025", "Median", ".975", "ESS", "Rhat")

      }

      overview_tabs[[i]] <- list(
        priors          = object$models[[i]]$priors,
        tab             = s.x,
        prior_prob      = prior_prob[i],
        marg_lik        = marg_lik[i],
        posterior_prob  = posterior_prob[i],
        BF              = BF[i],
        add_info        = list(
          BF_type      = paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")")
        )
      )
    }


    res <- list(
      call     = object$call,
      overview = overview_tabs,
      add_info = list(
        n_models         = length(object$models),
        digits_estimates = digits_estimates,
        digits_BF        = digits_BF,
        exp              = exp,
        BF_type          = paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")"),
        type             = "individual"
      )

    )

  }

  class(res) <- "summary.RoBTT"
  return(res)
}

#' @title Prints summary object for RoBTT method
#'
#' @param x a summary of a RoBTT object
#' @param ... additional arguments
#' @method print.summary RoBTT
#' @export print.summary.RoBTT
#' @rawNamespace S3method(print, summary.RoBTT)
#' @seealso [RoBTT()]
print.summary.RoBTT <- function(x, ...){

  # format the output before printing
  if(x$add_info$type == "ensemble"){

    overview <- x$overview
    overview$Models <- paste0(overview$Models,"/",  x$add_info$n_models - x$add_info$failed)
    overview[,2:3]  <- format(round(overview[,2:3], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
    overview[,4]    <- format(round(overview[,4],   x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)

    # round the results (a loop is required to deal with NAs)
    parameters_averaged <- x$parameters$averaged
    for(i in 1:ncol(parameters_averaged)){
      if(i == ncol(parameters_averaged)){
        parameters_averaged[,i] <- format(round(parameters_averaged[,i], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      }else{
        parameters_averaged[,i] <- format(round(parameters_averaged[,i], x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)
      }
    }
    estimates_averaged <- x$estimates$averaged
    for(i in 1:ncol(estimates_averaged)){
      if(i == ncol(estimates_averaged)){
        estimates_averaged[,i] <- format(round(estimates_averaged[,i], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      }else{
        estimates_averaged[,i] <- format(round(estimates_averaged[,i], x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)
      }
    }
    
    if(!is.null(x$parameters$conditional)){
      parameters_conditional <- x$parameters$conditional
      for(i in 1:ncol(parameters_conditional)){
        parameters_conditional[,i] <- format(round(parameters_conditional[,i], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      }
    }
    if(!is.null(x$estimates$conditional)){
      estimates_conditional <- x$estimates$conditional
      for(i in 1:ncol(estimates_conditional)){
        estimates_conditional[,i] <- format(round(estimates_conditional[,i], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      }
    }

  }else if(x$add_info$type == "models"){

    overview <- x$overview
    for(cn in c("Prior prob.",  "Post. prob.", "log(MargLik)")){
      overview[,cn]  <- format(round(overview[,cn], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
    }
    overview[,ncol(overview)]    <- format(round(overview[,ncol(overview)],   x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)

    if(!is.null(x$diagnostics)){
      diagnostics     <- x$diagnostics
      diagnostics[,"max(Rhat)"] <- format(round(diagnostics[,"max(Rhat)"], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      diagnostics[,"min(ESS)"]  <- format(round(diagnostics[,"min(ESS)"],  x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)
    }
  }else if(x$add_info$type == "individual"){

    overview <- x$overview

    for(i in 1:length(overview)){

      temp_main_info_names <- c("Model:", "Prior prob.:", "log(MargLik):", "Post. prob.:", paste0(x$add_info$BF_type,":"))
      temp_main_info       <- c(i,
                                format(round(overview[[i]]$prior_prob, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates),
                                format(round(overview[[i]]$marg_lik, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates),
                                format(round(overview[[i]]$posterior_prob, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates),
                                format(round(overview[[i]]$BF, x$add_info$digits_BF), nsmall = x$add_info$digits_BF))


      temp_main_priors_names <- "Parameter"
      temp_main_priors       <- "Prior Distributions"
      for(j in seq_along(overview[[i]]$priors)){
        temp_main_priors_names <- c(temp_main_priors_names, names(overview[[i]]$priors)[j])
        temp_main_priors       <- c(temp_main_priors, print(overview[[i]]$priors[[j]], silent = TRUE))
      }

      if(length(temp_main_info) > length(temp_main_priors)){
        temp_main_priors       <- c(temp_main_priors,       rep("", length(temp_main_info) - length(temp_main_priors)))
        temp_main_priors_names <- c(temp_main_priors_names, rep("", length(temp_main_info) - length(temp_main_priors_names)))
      }else if(length(temp_main_priors) > length(temp_main_info)){
        temp_main_info       <- c(temp_main_priors,       rep("", length(temp_main_priors) - length(temp_main_info)))
        temp_main_info_names <- c(temp_main_priors_names, rep("", length(temp_main_priors) - length(temp_main_info_names)))
      }

      temp_main <- data.frame(cbind(
        temp_main_info_names, temp_main_info,
        rep("", length(temp_main_info)), rep("", length(temp_main_info)),
        temp_main_priors_names, temp_main_priors))
      colnames(temp_main) <- c("", " ", "  ", "   ", "    ", "     ")
      overview[[i]]$tab_main <- temp_main

      if(!is.null(overview[[i]]$tab)){
        temp_tab <- overview[[i]]$tab
        overview[[i]]$tab[, -6] <- format(round(temp_tab[,-6], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
        overview[[i]]$tab[,  6] <- round(temp_tab[,6])
      }

      overview[[i]]$prior_prob     <- format(round(overview[[i]]$prior_prob, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      overview[[i]]$marg_lik       <- format(round(overview[[i]]$marg_lik, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      overview[[i]]$posterior_prob <- format(round(overview[[i]]$posterior_prob, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      overview[[i]]$BF             <- format(round(overview[[i]]$BF, x$add_info$digits_BF), nsmall = x$add_info$digits_BF)
    }
  }


  cat("Call:\n")
  print(x$call)
  cat("\n")


  if(x$add_info$type == "ensemble"){

    cat("Robust Bayesian T-Test\n")
    print(overview, quote = FALSE, right = TRUE)
    cat("\n")

    cat("Model-averaged parameter estimates\n")
    print(parameters_averaged, quote = FALSE, right = TRUE)
    
    if(!is.null(estimates_averaged)){
      cat("\n")
      cat("Model-averaged group estimates\n")
      print(estimates_averaged, quote = FALSE, right = TRUE)
    }

    if(!is.null(x$parameters$conditional)){
      cat("\n")
      cat("Conditional parameter estimates\n")
      print(parameters_conditional, quote = FALSE, right = TRUE)
    }
    
    if(!is.null(x$estimates$conditional)){
      cat("\n")
      cat("Conditional group estimates\n")
      print(estimates_conditional, quote = FALSE, right = TRUE)
    }

    if(x$add_info$failed != 0)cat(paste0("\033[0;31m",x$add_info$failed, ifelse(x$add_info$failed == 1, " model", " models"), " failed to converge and ",ifelse(x$add_info$failed == 1, "was", "were")," omited from the summary.\033[0m\n"))
    
  }else if(x$add_info$type == "models"){

    cat("Robust Bayesian T-Test\n")
    print(overview, quote = FALSE, right = TRUE)

    if(!is.null(x$diagnostics)){
      cat("\n")
      cat("Models diagnostics overview\n")
      print(diagnostics, quote = FALSE, right = TRUE)
    }

  }else if(x$add_info$type == "individual"){

    cat("Individual Models Summary\n\n")

    for(i in 1:length(overview)){

      cat("Model Overview:\n")
      print(overview[[i]]$tab_main, quote = FALSE, right = TRUE, row.names = FALSE)

      cat("\nModel Coefficients:\n")
      print(overview[[i]]$tab, quote = FALSE, right = TRUE)

      if(i != length(overview))cat("\n\n")
    }
  }

}

#' @title Reports whether x is a RoBTT object
#'
#' @param x an object to test
#' @export is.RoBTT
is.RoBTT            <- function(x){
  inherits(x, "RoBTT")
}

