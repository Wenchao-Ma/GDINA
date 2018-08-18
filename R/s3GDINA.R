#' @include GDINA.R


#' @export
#' @describeIn GDINA calculate log-likelihood
logLik.GDINA <- function (object, ...){

  logL <- extract(object, what = "logLik")

  attr(logL, "df") <- extract(object, what = "npar")

  attr(logL, "nobs") <- extract(object, what = "nobs")

  class(logL) <- "logLik"

  return(logL)

}

#' @export
#' @describeIn GDINA calculate deviance
deviance.GDINA <- function(object,...){
  extract(object, what = "deviance")
}
#' @export
#' @describeIn GDINA calculate number of observations
nobs.GDINA <- function(object,...){
  extract(object, what = "nobs")
}

#' @export
#' @describeIn GDINA calculate covariance-matrix for delta parameters
vcov.GDINA <- function(object,...){
  dn <- lapply(coef(object,"delta"),names)
  nm <- names(dn)
  v.names <- NULL
  for(j in seq_len(length(dn))) v.names <- c(v.names,paste(nm[j],dn[[j]]))
  opg <- OPG_d(object,SE.type = 2)$cov
  rownames(opg) <- colnames(opg) <- v.names
  opg
  }

#' @title Calculate the number of parameters
#' @description Calculate the number of parameters for GDINA estimates.
#' Returned the total number of parameters, the number of item parameters and
#' the number parameters of joint attribute distribution.
#' @param object GDINA object
#' @param ... additional arguments
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#'
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' npar(fit)
#'}
#' @export
#'
npar <- function(object,...){
  UseMethod("npar")
}
#' @export
#' @describeIn GDINA calculate the number of parameters
npar.GDINA <- function(object,...){
  out <- list(npar=extract(object, what = "npar"),
    npar.item=extract(object, what = "npar.item"), # free param
    npar.att=extract(object, what = "npar.att"),
    npar.total.item=extract(object, what = "npar.item") + extract(object, what = "npar.fixeditem"),
    npar.fix.item=extract(object, what = "npar.fixeditem"))
  names(out) <- c("No. of parameters","No. of estimated item parameters","No. of population parameters",
                  "No. of total item parameters","No. of fixed item parameters")
  class(out) <- "npar.GDINA"
  return(out)
}


#' @title Extract log-likelihood for each individual
#' @description Extract individual log-likelihood.
#' @param object GDINA object
#' @param ... additional arguments
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#'
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' iL <- indlogLik(fit)
#' iL[1:6,]
#'}
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
#' @description Extract individual log posterior.
#' @param object GDINA object
#' @param ... additional arguments
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' iP <- indlogPost(fit)
#' iP[1:6,]
#'}
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

