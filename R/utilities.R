#' @title Options for the 'RoBTT' package
#'
#' @description A placeholder object and functions for the 'RoBTT' package.
#' (adapted from the runjags R package).
#'
#' @param name the name of the option to get the current value of - for a list of
#' available options, see details below.
#' @param ... named option(s) to change - for a list of available options, see
#' details below.
#'
#' @return The current value of all available 'RoBTT' options (after applying any
#' changes specified) is returned invisibly as a named list.
#'
#' @export RoBTT.options
#' @export RoBTT.get_option
#' @name RoBTT_options
#' @aliases RoBTT_options RoBTT.options RoBTT.get_option
NULL


#' @rdname RoBTT_options
RoBTT.options    <- function(...){
  
  opts <- list(...)
  
  for(i in seq_along(opts)){
    
    if(!names(opts)[i] %in% names(RoBTT.private))
      stop(paste("Unmatched or ambiguous option '", names(opts)[i], "'", sep=""))
    
    assign(names(opts)[i], opts[[i]] , envir = RoBTT.private)
  }
  
  return(invisible(RoBTT.private$options))
}

#' @rdname RoBTT_options
RoBTT.get_option <- function(name){
  
  if(length(name)!=1)
    stop("Only 1 option can be retrieved at a time")
  
  if(!name %in% names(RoBTT.private))
    stop(paste("Unmatched or ambiguous option '", name, "'", sep=""))
  
  # Use eval as some defaults are put in using 'expression' to avoid evaluating at load time:
  return(eval(RoBTT.private[[name]]))
}



# adapted from the runjags package version 2.2.0
RoBTT.private <- new.env()
# Use 'expression' for functions to avoid having to evaluate before the package is fully loaded:
assign("defaultoptions",  list(
  jagspath = expression(runjags::findjags()),
  envir    = RoBTT.private))

assign("options",         RoBTT.private$defaultoptions,   envir = RoBTT.private)
assign("RoBTT_version",   "notset",                       envir = RoBTT.private)
assign("max_cores",       parallel::detectCores(logical = TRUE) - 1,  envir = RoBTT.private)

# check proper BayesTools package version
.check_BayesTools <- function(){
  
  RoBTT.version      <- try(utils::packageVersion("RoBTT"))
  BayesTools.version <- try(utils::packageVersion("BayesTools"))
  
  if(inherits(RoBTT.version, "try-error") | inherits(BayesTools.version, "try-error")){
    return(invisible(FALSE))
  }
  
  if(is.null(RoBTT.version) | is.null(BayesTools.version)){
    return(invisible(FALSE))
  }
  
  BayesTools_required <- switch(
    paste0(RoBTT.version, collapse = "."),
    "1.0.0" = c("0.2.12", "999.999.999"),
    "1.0.1" = c("0.2.12", "999.999.999"),
    "1.0.2" = c("0.2.12", "999.999.999"),
    "1.0.3" = c("0.2.12", "999.999.999"),
    "1.1.0" = c("0.2.14", "999.999.999"),
    stop("New RoBTT version needs to be defined in '.check_BayesTools' function!")
  )
  
  min_OK <- sum(as.numeric(strsplit(BayesTools_required[1], ".", fixed = TRUE)[[1]]) * c(1e9, 1e6, 1e3)) <=
    sum(unlist(BayesTools.version) * c(1e9, 1e6, 1e3))
  max_OK <- sum(as.numeric(strsplit(BayesTools_required[2], ".", fixed = TRUE)[[1]]) * c(1e9, 1e6, 1e3)) >=
    sum(unlist(BayesTools.version) * c(1e9, 1e6, 1e3))
  
  if(min_OK && max_OK){
    return(invisible(TRUE))
  }else{
    warning(sprintf(
      "RoBTT version %1$s requires BayesTools version higher or equal %2$s and lower or equal %3$s.",
      paste0(RoBTT.version, collapse = "."),
      BayesTools_required[1],
      BayesTools_required[2]
    ), call.= FALSE)
    return(invisible(FALSE))
  }
}