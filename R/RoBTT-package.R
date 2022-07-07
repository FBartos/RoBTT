##' RoBTT: Robust Bayesian t-test
##'
##' RoBTT: Bayesian model-averaged t-test extends the Bayesian t-test by 
##' incorporating inference about heterogeneity of variances and outliers.
##'
##'
##' @name RoBTT-package
##' @author František Bartoš \email{f.bartos96@@gmail.com}
##' @keywords package
##' @aliases RoBTT-package RoBTT_package RoBTT.package
##' @docType package
##' @section
##' User guide: See \insertCite{maier2022bayesian;textual}{RoBTT} for 
##' details regarding the RoBTT methodology.
##' 
##' More details regarding customization of the Bayesian model-averaged t-test 
##' are provided in 
##' \href{../doc/Introduction_to_RoBTT.html}{\bold{Introduction to RoBTT}} 
##' vignette. Please, use the "Issues" section in the GitHub repository to 
##' ask any further questions.
##'
##' @references \insertAllCited{}
##' @importFrom BayesTools is.prior is.prior.none is.prior.point is.prior.simple
##' @rawNamespace import(Rcpp)
##' @rawNamespace useDynLib(RoBTT, .registration = TRUE)
"_PACKAGE"
