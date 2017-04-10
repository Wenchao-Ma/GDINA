#' @include GDINA.R

#'@export
#'@describeIn GDINA calculate AIC
AIC.GDINA <- function(object,...){
  aic <- extract(object, what = "AIC")
  class(aic) <- "AIC.GDINA"
  aic
}

#' @export
#' @describeIn GDINA calculate BIC
BIC.GDINA <- function(object,...){
  bic <- extract(object, what = "BIC")
  class(bic) <- "BIC.GDINA"
  bic
}

#' @export
#' @describeIn GDINA calculate log-likelihood
logLik.GDINA <- function(object,...){
  logL <- extract(object, what = "logLik")
  class(logL) <- "logLik.GDINA"
  logL
}
#' @export
#' @describeIn GDINA calculate deviance
deviance.GDINA <- function(object,...){
  dev <- extract(object, what = "deviance")
  class(dev) <- "deviance.GDINA"
  dev
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
  out <- list(npar=extract(object, what = "npar"),
    npar.item=extract(object, what = "npar.item"),
    npar.att=extract(object, what = "npar.att"))
  names(out) <- c("No. of parameters","No. of item parameters","No. of population parameters")
  class(out) <- "npar.GDINA"
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
