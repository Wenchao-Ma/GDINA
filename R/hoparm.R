#'@include GDINA-package.R GDINA.R
#'@title extract higher-order parameters
#'
#' @description
#' Function to extract higher-order parameters when a higher-order model is fitted.
#'
#' @param object estimated GDINA object returned from \code{\link{GDINA}}
#' @param withSE estimate standard errors for lambda parameters or not?
#' @param digits how many decimal places for the ouput?
#' @param ... additional arguments
#' @return a list with element \code{theta} for higher-order incidental (ability) parameters
#' and \code{lambda} for higher-order structural parameters.
#'
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#'
#'@export
hoparm <- function(object, withSE = FALSE,digits = 4, ...) {
  UseMethod("hoparm")
}


#' @title NULL
#' @description To extract higher-order parameters, use method \code{\link{hoparm}}.
#' @describeIn GDINA extract higher-order parameters
#' @aliases itemparm.GDINA
#' @export
hoparm.GDINA <- function(object, withSE = FALSE,digits = 4, ...){
  if(!internalextract(object,what = "higher.order"))stop("Set higher.order = TRUE to estimate a higher-order model.",call. = FALSE)
  K <- extract.GDINA(object,what = "natt")
  Q <- extract.GDINA(object,what = "Q")
  pattern <- alpha(K,T,Q)
  theta <- lambda <- NULL
  quad <- seq(-4,4,by=0.1)
  lk <- HO.loglik(extract.GDINA(object,what = "higher.order.struc.parm")[,1],
                  extract.GDINA(object,what = "higher.order.struc.parm")[,2],
                  theta=quad,X=pattern) #nnode x 2^K
  theta <- round(t(apply(exp(extract.GDINA(object,what = "logposterior.i")),1,function(x){
    est <- sum(quad*colSums(exp(t(lk))*x)*dnorm(quad))/
      sum(colSums(exp(t(lk))*x)*dnorm(quad))
    se <- sqrt(sum((quad-est)^2*colSums(exp(t(lk))*x)*dnorm(quad))/
                 sum(colSums(exp(t(lk))*x)*dnorm(quad)))
    return(c(est,se))

  })),digits)
  colnames(theta) <- c("ability","S.E.")


    if(withSE) {
      ho <- internalextract(object,what = "higher.order.struc.parm")
      K=internalextract(object,what = "natt")
      HO.se <- NULL
      if (internalextract(object,what = "higher.order.model")=="2PL") {
        inv.info <- (solve((-1)*numDeriv::hessian(func=HO.SE.2,x=c(ho[,1],ho[,2]),
                                                  Xloglik=internalextract(object,what = "loglikelihood.i"),
                                                  K=K,nnodes=19,
                                                  N=internalextract(object,what = "nobs"))))
        HO.se <- matrix(sqrt(diag(inv.info)),ncol = 2)
        colnames(HO.se) <- c("slope.se","intercept.se")
      }else if (internalextract(object,what = "higher.order.model")=="1PL") {
        inv.info <- (solve((-1)*numDeriv::hessian(func=HO.SE.1,x=c(ho[,1],ho[,2]),
                                                  Xloglik=internalextract(object,what = "loglikelihood.i"),
                                                  K=K,nnodes=19,
                                                  N=internalextract(object,what = "nobs"))))
        HO.se <- sqrt(diag(inv.info))
        HO.se <- cbind(rep(HO.se[1],K),HO.se[2:(K+1)])
        colnames(HO.se) <- c("slope.se","intercept.se")
      }else if (internalextract(object,what = "higher.order.model")=="Rasch") {
        inv.info <- (solve((-1)*numDeriv::hessian(func=HO.SE.R,x=c(ho[,2]),
                                                  Xloglik=internalextract(object,what = "loglikelihood.i"),
                                                  K=K,nnodes=19,
                                                  N=internalextract(object,what = "nobs"))))
        HO.se <- cbind(NA,sqrt(diag(inv.info)))
        colnames(HO.se) <- c("slope.se","intercept.se")
      }
      lambda <- round(cbind(ho,HO.se),digits)

    }else{
      lambda <- round(internalextract(object,what = "higher.order.struc.parm"),digits)
    }


  return(list(theta=theta,lambda=lambda))

}
