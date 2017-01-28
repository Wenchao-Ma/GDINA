#'@include GDINA-package.R GDINA.R
#'@title calculate lower-order incidental (person) parameters
#'
#'
#' @description
#' Function to calculate various person attribute parameters, including \code{"EAP"},
#' \code{"MAP"}, and \code{"MLE"}, for EAP, MAP and MLE estimates of
#' attribute patterns, \code{"mp"} for marginal mastery probabilities.
#' See \code{\link{GDINA}} for examples. To estimate higher-order person parameters,
#' see \code{\link{hoparm}}.
#'
#'
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @param object estimated GDINA object returned from \code{\link{GDINA}}
#' @param what what to extract; It can be \code{"EAP"},
#' \code{"MAP"}, and \code{"MLE"}, for EAP, MAP and MLE estimates of
#' attribute patterns, and \code{"mp"} for marginal mastery probabilities.
#' @param digits number of decimal places.
#' @param ... additional arguments
#'
#' @references
#' Huebner, A., \& Wang, C. (2011). A note on comparing examinee classification methods for cognitive diagnosis models. \emph{Educational and Psychological Measurement, 71}, 407-419.
#'
#'
#'@export
personparm <- function (object, what=c("EAP","MAP","MLE", "mp"),digits = 4,...) {
  UseMethod("personparm")
}

#' @title NULL
#' @description To calculate lower-order incidental (person) parameters
#' use method \code{\link{personparm}}. To extract other components returned, use \code{\link{extract}}.
#' To plot item/category response function, use \code{\link{plotIRF}}. To
#' check whether monotonicity is violated, use \code{\link{monocheck}}. To conduct anaysis in graphical user interface,
#' use \code{\link{startGDINA}}.
#' @describeIn GDINA calculate person attribute patterns and higher-order ability
#' @aliases personparm.GDINA
#' @export
personparm.GDINA <- function(object,
                             what=c("EAP","MAP","MLE", "mp"),digits = 4,...){
  what <- match.arg(what)
  # The number of attributes
  K <- internalextract(object,what = "natt")
  # Q-matrix
  Q <- internalextract(object,what = "Q")
  pattern <- alpha(K,T,Q)
  out <- NULL
      # dichotomous attributes
      switch(what,
             EAP={
        if (max(Q) == 1)
        {
          # dichotomous attributes
          out <- 1*((exp(internalextract(object,what = "logposterior.i")) %*% pattern) > 0.5000)
        }else{
          # polytomous attributes
          attnum <- list()
          for (k in 1:K)
          {
            tmp <- NULL
            for (kj in 0:max(Q[, k]))
            {
              tmp <- cbind(tmp, apply(exp(internalextract(object,what = "logposterior.i"))[, which(pattern[, k] == kj),drop = FALSE], 1, sum))
            }

            attnum[[k]] <- tmp
            out <- cbind(out, max.col(tmp))
          }

          out <- out - 1

        }
        colnames(out)[1:K] <- paste("A",1:K,sep = "")
      },
      MLE={
        if(max(Q) > 1) stop("Please use EAP for polytomous attributes.", call. = FALSE)
        out <- data.frame(MLE=pattern[max.col(internalextract(object,what = "loglikelihood.i")),],
                          multimodes=as.logical(rowSums(internalextract(object,what = "loglikelihood.i")==apply(internalextract(object,what = "loglikelihood.i"),1,max))-1))
        colnames(out)[1:K] <- paste("A",1:K,sep = "")
      },
MAP={
        if(max(Q) > 1) stop("Please use EAP for polytomous attributes.", call. = FALSE)
        out <- data.frame(MAP=pattern[max.col(internalextract(object,what = "logposterior.i")),],
                          multimodes=as.logical(rowSums(internalextract(object,what = "logposterior.i")==apply(internalextract(object,what = "logposterior.i"),1,max))-1))
        colnames(out)[1:K] <- paste("A",1:K,sep = "")
      },
mp={
        if(max(Q) > 1) stop("Not available for polytomous attributes.", call. = FALSE)
          out <- round(exp(internalextract(object,what = "logposterior.i")) %*% pattern,
                       digits)
          colnames(out)[1:K] <- paste("A",1:K,sep = "")
        })

return(out)

}
