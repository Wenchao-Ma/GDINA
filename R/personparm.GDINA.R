#'@include GDINA-package.R GDINA.R
#'@title calculate person (incidental) parameters
#'
#'
#' @description
#' Function to calculate various person attribute parameters, including \code{"EAP"},
#' \code{"MAP"}, and \code{"MLE"}, for EAP, MAP and MLE estimates of
#' attribute patterns (see Huebner & Wang, 2011), \code{"mp"} for marginal mastery probabilities, and \code{"HO"}
#' for higher-order ability estimates if a higher-order model is fitted.
#' See \code{\link{GDINA}} for examples.
#'
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @param object estimated GDINA object returned from \code{\link{GDINA}}
#' @param what what to extract; It can be \code{"EAP"},
#' \code{"MAP"}, and \code{"MLE"}, for EAP, MAP and MLE estimates of
#' attribute patterns, and \code{"mp"} for marginal mastery probabilities, and \code{"HO"}
#' for higher-order ability estimates if a higher-order model is fitted.
#' @param digits number of decimal places.
#' @param ... additional arguments
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' # EAP
#' head(personparm(fit))
#' # MAP
#' head(personparm(fit, what = "MAP"))
#'}
#'
#' @references
#' Huebner, A., & Wang, C. (2011). A note on comparing examinee classification methods for cognitive diagnosis models. \emph{Educational and Psychological Measurement, 71}, 407-419.
#'
#'
#'@export
personparm <- function (object, what=c("EAP","MAP","MLE", "mp", "HO"),digits = 4,...) {
  UseMethod("personparm")
}

#' @title NULL
#' @description To calculate lower-order incidental (person) parameters
#' use method \code{\link{personparm}}. To extract other components returned, use \code{\link{extract}}.
#' To plot item/category response function, use \code{\link{plot}}. To
#' check whether monotonicity is violated, use \code{\link{monocheck}}. To conduct anaysis in graphical user interface,
#' use \code{\link{startGDINA}}.
#' @param object GDINA object for various S3 methods
#' @param what argument for various S3 methods; For calculating structural parameters using \code{\link{coef}},
#' \code{what} can be
#'  \itemize{
#'    \item \code{itemprob} - item success probabilities of each reduced attribute pattern.
#'    \item \code{catprob} - category success probabilities of each reduced attribute pattern; the same as \code{itemprob} for dichtomous response data.
#'    \item \code{LCprob} - item success probabilities of each attribute pattern.
#'    \item \code{gs} - guessing and slip parameters of each item/category.
#'    \item \code{delta} - delta parameters of each item/category, see G-DINA formula in details.
#'    \item \code{rrum} - RRUM parameters when items are estimated using \code{RRUM}.
#'    \item \code{lambda} - structural parameters for joint attribute distribution.
#'    }
#'    For calculating incidental parameters using \code{\link{personparm}},
#' \code{what} can be
#'  \itemize{
#'    \item \code{EAP} - EAP estimates of attribute pattern.
#'    \item \code{MAP} - MAP estimates of attribute pattern.
#'    \item \code{MLE} - MLE estimates of attribute pattern.
#'    \item \code{mp} -  marginal mastery probabilities.
#'    \item \code{HO} -  EAP estimates of higher-order ability if a higher-order is fitted.
#'    }
#' @param withSE argument for method \code{\link{coef}}; estimate standard errors or not?
#' @param SE.type type of standard errors. For now, SEs are calculated based on outper-product of gradient.
#'   It can be \code{1} based on item-wise information, \code{2} based on incomplete information and \code{3}
#'   based on complete information.
#' @param digits How many decimal places in each number? The default is 4.
#' @describeIn GDINA calculate person attribute patterns and higher-order ability
#' @aliases personparm.GDINA
#' @export
personparm.GDINA <- function(object,
                             what=c("EAP","MAP","MLE", "mp", "HO"),digits = 4,...){
  what <- match.arg(what)
  # The number of attributes
  K <- extract(object,what = "natt")
  # Q-matrix
  Q <- extract(object,what = "Q")
  pattern <- attributepattern(Q=Q)
  out <- NULL
      # dichotomous attributes
      switch(what,
             EAP={
        if (max(Q) == 1)
        {
          # dichotomous attributes
          out <- 1*((exp(extract(object,what = "logposterior.i")) %*% pattern) > 0.5000)
        }else{
          # polytomous attributes
          attnum <- list()
          for (k in 1:K)
          {
            tmp <- NULL
            for (kj in 0:max(Q[, k]))
            {
              tmp <- cbind(tmp, apply(exp(extract(object,what = "logposterior.i"))[, which(pattern[, k] == kj),drop = FALSE], 1, sum))
            }

            attnum[[k]] <- tmp
            out <- cbind(out, max.col(tmp))
          }

          out <- out - 1

        }
        colnames(out)[1:K] <- paste("A",1:K,sep = "")
      },
      MLE={
        out <- data.frame(MLE=pattern[max.col(extract(object,what = "loglikelihood.i")),],
                          multimodes=as.logical(rowSums(extract(object,what = "loglikelihood.i")==apply(extract(object,what = "loglikelihood.i"),1,max))-1))
        colnames(out)[1:K] <- paste("A",1:K,sep = "")
      },
MAP={
        out <- data.frame(MAP=pattern[max.col(extract(object,what = "logposterior.i")),],
                          multimodes=as.logical(rowSums(extract(object,what = "logposterior.i")==apply(extract(object,what = "logposterior.i"),1,max))-1))
        colnames(out)[1:K] <- paste("A",1:K,sep = "")
      },
mp={
  post <- exp(extract(object,what = "logposterior.i"))
        if(max(Q) == 1){
          out <- round(post %*% pattern, digits)
          colnames(out)[1:K] <- paste("A",1:K,sep = "")
        }else{
          C <- sum(apply(pattern,2,max))
          out <- matrix(0,extract(object,"nobs"),C)
          nam <- vector("character",C)
          kkk <- 1
          for(k in 1:K){
            for(kk in 1:max(pattern[,k])){
              out[,kkk] <- rowSums(post[,which(pattern[,k]==kk),drop=FALSE])
              nam[kkk] <- paste0("A",k,"[Level",kk,"]",collapse = "")
              kkk <- kkk + 1
            }
          }
          colnames(out) <- nam
        }
         out
        },
HO={
  if(any(extract(object,what ="att.dist")!= "higher.order")) {
    stop("Set att.dist = 'higher.order' to estimate a higher-order model.",call. = FALSE)
  }


  quad <- extract(object,"higher.order")$QuadNodes
  w <- extract(object,"higher.order")$QuadWghts
  lambda <- NULL
  K <- extract(object,what = "natt")
  Q <- extract(object,what = "Q")
  N <- extract(object,"nobs")
  pattern <- attributepattern(K)
  Lx <- exp(extract(object,what = "logposterior.i"))
  g <- extract(object,"gr")

  post <- list()
  for(gg in 1:extract(object,"ngroup")){
    post[[gg]] = PostTheta(AlphaPattern = pattern, theta = quad, f_theta = w, a=extract(object,what = "struc.parm")[[gg]][,1],
                           b=extract(object,what = "struc.parm")[[gg]][,2]) # 2^K x nnodes P(theta_q|AlphaPattern)
  }

    theta <- se <- rep(0,N)

    for(i in 1:N){
      ptheta_i <- colSums(post[[g[i]]]*Lx[i,]) #vector of length nnodes
      theta[i] <- sum(quad[,g[i]]*ptheta_i)
      se[i] <- sqrt(sum((quad[,g[i]]-theta[i])^2*ptheta_i))
    }
    out <- round(data.frame(theta=theta,se=se),digits)
    colnames(out) <- c("EAP","SE")

})

return(out)

}
