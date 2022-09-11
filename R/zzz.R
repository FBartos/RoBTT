.onLoad <- function(libname, pkgname) {
  
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  
  for(m in modules){
    loadModule(m, what = TRUE)
  }
  
  RoBTT.private$RoBTT_version <- utils::packageDescription(pkgname, fields = 'Version')
  setopts <- mget('.RoBTT.options', envir=.GlobalEnv, ifnotfound = list(.RoBTT.options = NULL))[[1]]
  .check_BayesTools()
}
