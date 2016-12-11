#' @include GDINA.R

#'@export
#'@describeIn GDINA calculate AIC
AIC.GDINA <- function(object,...){
  extract.GDINA(object, what = "AIC")
}

#' @export
#' @describeIn GDINA calculate BIC
BIC.GDINA <- function(object,...){
  extract.GDINA(object, what = "BIC")
}

#' @export
#' @describeIn GDINA calculate log-likelihood
logLik.GDINA <- function(object,...){
  extract.GDINA(object, what = "logLik")
}
#' @export
#' @describeIn GDINA calculate deviance
deviance.GDINA <- function(object,...){
  extract.GDINA(object, what = "deviance")
}



#' @title Calculate the number of parameters
#' @description Calculate the number of parameters for GDINA estimates.
#' Returned the total number of parameters, the number of item parameters and
#' the number of population parameters.
#' See \code{\link{GDINA}} for examples.
#' @param object GDINA object
#' @param ... additional arguments
#' @export
#'
npar <- function(object,...){
  UseMethod("npar")
}
#' @export
#' @describeIn GDINA calculate the number of parameters
npar.GDINA <- function(object,...){
  out <- c(extract.GDINA(object, what = "npar"),
    extract.GDINA(object, what = "npar.item"),
    extract.GDINA(object, what = "npar.att"))
  names(out) <- c("# of parameters","# of item parameters","# of population parameters")
  return(out)
}


#' @title Extract log-likelihood for each individual
#' @description Extract individual log-likelihood
#' See \code{\link{GDINA}} for examples.
#' @param object GDINA object
#' @param ... additional arguments
#' @export
#'
indlogLik <- function(object,...){
  UseMethod("indlogLik")
}
#' @export
#' @describeIn GDINA extract log-likelihood for each individual
indlogLik.GDINA <- function(object,...){
  return(extract(object,"loglikelihood.i"))
}

#' @title Extract log posterior for each individual
#' @description Extract individual log posterior
#' See \code{\link{GDINA}} for examples.
#' @param object GDINA object
#' @param ... additional arguments
#' @export
#'
indlogPost <- function(object,...){
  UseMethod("indlogPost")
}
#' @export
#' @describeIn GDINA extract log posterior for each individual
indlogPost.GDINA <- function(object,...){
  return(extract(object,"logposterior.i"))
}
