.onLoad <- function(libname, pkgname) {
  RoBTT.private$RoBTT_version <- utils::packageDescription(pkgname, fields = 'Version')
  setopts <- mget('.RoBTT.options', envir=.GlobalEnv, ifnotfound = list(.RoBTT.options = NULL))[[1]]
  .check_BayesTools()
}
