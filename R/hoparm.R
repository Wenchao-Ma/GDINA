#'@include GDINA-package.R GDINA.R
#'@title extract higher-order parameters
#'
#' @description
#' Function to extract higher-order parameters when a higher-order model is fitted.
#'
#' @param object estimated GDINA object returned from \code{\link{GDINA}}
#' @param theta.est logical; Estimating higher-order person ability or not? The default is \code{FALSE}.
#' @param withSE estimate standard errors for lambda parameters or not?
#' @param digits how many decimal places for the ouput?
#' @param ... additional arguments
#' @return a list with element \code{theta} for higher-order incidental (ability) parameters
#' and \code{lambda} for higher-order structural parameters.
#'
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#'
#'@export
hoparm <- function(object, withSE = FALSE,theta.est = FALSE, digits = 4, ...) {
  UseMethod("hoparm")
}


#' @title NULL
#' @description To extract higher-order parameters, use method \code{\link{hoparm}}.
#' @param theta.est logical; Estimating higher-order person ability or not? The default is \code{FALSE}.
#' @describeIn GDINA extract higher-order parameters
#' @aliases itemparm.GDINA
#' @export
hoparm.GDINA <- function(object, withSE = FALSE,theta.est = FALSE, digits = 4, ...){
  if(!class(object)=="GDINA") stop("object must be a GDINA estimate.",call. = FALSE)
  if(all(extract(object,what ="att.dist")!= "higher.order")) {
    stop("Set att.dist = 'higher.order' to estimate a higher-order model.",call. = FALSE)
  }else{
    if(extract(object,"ngroup")>1) {
      withSE <- theta.est <- FALSE
    warning("SE for higher-order structural parameters cannot be estimated for multiple groups",call. = FALSE)}
  }

  theta <- NULL
  quad <- seq(-4,4,length.out = extract(object,"higher.order")$nquad)

  if(theta.est){
    K <- extract(object,what = "natt")
    Q <- extract(object,what = "Q")
    pattern <- alpha(K,T,Q)
    theta <- lambda <- NULL
    lk <- HO.loglik(extract(object,what = "higher.order.struc.parm")[,1],
                    extract(object,what = "higher.order.struc.parm")[,2],
                    theta=quad,X=pattern) #nnode x 2^K
    theta <- round(t(apply(exp(extract(object,what = "logposterior.i")),1,function(x){
      est <- sum(quad*colSums(exp(lk)*x)*dnorm(quad))/
        sum(colSums(exp(lk)*x)*dnorm(quad))
      se <- sqrt(sum((quad-est)^2*colSums(exp(lk)*x)*dnorm(quad))/
                   sum(colSums(exp(lk)*x)*dnorm(quad)))
      return(c(est,se))

    })),digits)
    colnames(theta) <- c("ability","S.E.")

  }


    if(withSE) {
      ho <- extract(object,what = "higher.order.struc.parm")
      K=extract(object,what = "natt")
      HO.se <- NULL
      if (extract(object,what = "higher.order.model")=="2PL") {
        inv.info <- (solve((-1)*numDeriv::hessian(func=HO.SE.2,x=c(ho[,1],ho[,2]),
                                                  Xloglik=extract(object,what = "loglikelihood.i"),
                                                  K=K,nnodes=length(quad),
                                                  N=extract(object,what = "nobs"))))
        HO.se <- matrix(sqrt(diag(inv.info)),ncol = 2)
        colnames(HO.se) <- c("slope.se","intercept.se")
      }else if (extract(object,what = "higher.order.model")=="1PL") {
        inv.info <- (solve((-1)*numDeriv::hessian(func=HO.SE.1,x=c(ho[,1],ho[,2]),
                                                  Xloglik=extract(object,what = "loglikelihood.i"),
                                                  K=K,nnodes=length(quad),
                                                  N=extract(object,what = "nobs"))))
        HO.se <- sqrt(diag(inv.info))
        HO.se <- cbind(rep(HO.se[1],K),HO.se[2:(K+1)])
        colnames(HO.se) <- c("slope.se","intercept.se")
      }else if (extract(object,what = "higher.order.model")=="Rasch") {
        inv.info <- (solve((-1)*numDeriv::hessian(func=HO.SE.R,x=c(ho[,2]),
                                                  Xloglik=extract(object,what = "loglikelihood.i"),
                                                  K=K,nnodes=length(quad),
                                                  N=extract(object,what = "nobs"))))
        HO.se <- cbind(NA,sqrt(diag(inv.info)))
        colnames(HO.se) <- c("slope.se","intercept.se")
      }
      lambda <- round(cbind(ho,HO.se),digits)

    }else{
      if(extract(object,"ngroup")==1){
        lambda <- round(extract(object,what = "higher.order.struc.parm"),digits)
      }else{
        lambda <- lapply(extract(object,what = "higher.order.struc.parm"),function(x) if(!is.null(x)) round(x,digits = digits))
      }

    }


  return(list(theta=theta,lambda=lambda))

}
