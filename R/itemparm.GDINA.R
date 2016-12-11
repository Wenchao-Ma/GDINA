#'@include GDINA-package.R GDINA.R
#'@title extract lower-order structural (item) parameters
#'
#' @description
#' Function to extract various item parameters, including \code{"itemprob"} for
#' item success probabilities of each reduced attribute pattern, \code{"catprob"} for
#' category success probabilities of each reduced attribute pattern, \code{"LCprob"} for
#' item success probabilities of each attribute pattern,\code{"gs"} for guessing and slip parameters,
#' \code{"delta"} for delta parameters, \code{"rrum"} for RRUM parameters when items
#' are estimated using RRUM and \code{"higher.order"} for higher-order IRT
#' parameters. Standard errors can be estimated if \code{withSE = TRUE}.
#' See \code{\link{GDINA}} for examples.
#'
#' @param object estimated GDINA object returned from \code{\link{GDINA}}
#' @param what what to show; It can be \code{"itemprob"} for
#' item success probabilities of each reduced attribute pattern, \code{"catprob"} for
#' category success probabilities of each reduced attribute pattern, \code{"LCprob"} for
#' item success probabilities of each attribute pattern, \code{"gs"} for guessing and slip parameters,
#' \code{"delta"} for delta parameters, \code{"rrum"} for RRUM parameters when items
#' are estimated using RRUM and \code{"higher.order"} for higher-order IRT
#' parameters. The default is \code{"itemprob"}.
#' @param withSE show standard errors or not?
#' @param SE.type type of standard errors.
#' @param digits how many decimal places for the ouput?
#' @param ... additional arguments
#'
#'
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#'
#'@export
itemparm <- function(object,
                      what = c("catprob","gs","delta","rrum","itemprob","LCprob"),
                      withSE = FALSE, SE.type = 2,digits = 4, ...) {
  UseMethod("itemparm")
}


#' @title NULL
#' @description To extract lower-order structural (item) parameters, use S3 method \code{\link{itemparm}}.
#' @param object estimated GDINA object for various S3 methods
#' @param what argument for various S3 methods
#' @param withSE argument for S3 method \code{\link{itemparm}}; show standard errors or not?
#' @param SE.type type of standard errors.
#' @param digits How many decimal places in each number? The default is 4.
#' @param ... additional arguments
#' @describeIn GDINA extract various item parameters
#' @aliases itemparm.GDINA
#' @export
itemparm.GDINA <- function(object,
                           what = c("catprob","itemprob","LCprob","gs","delta","rrum","higher.order"),
                           withSE = FALSE, SE.type = 2,digits = 4, ...){
  what <- match.arg(what)
  if(tolower(what)=="catprob"){
    if(withSE){
      out <- mapply(rbind,internalextract(object,what = "catprob.parm"),
                    internalextract(object,what = "catprob.se", type = SE.type),SIMPLIFY = F)
      out <- lapply(out,function(x){rownames(x) <- c("Est.","S.E.");round(x,digits)})
    }else{
      out <- lapply(internalextract(object,what = "catprob.parm"),round,digits)
    }
  }else if(tolower(what)=="itemprob"){
    if(withSE&!internalextract(object,what = "sequential")){
        out <- mapply(rbind,internalextract(object,what = "itemprob.parm"),
                      internalextract(object,what = "itemprob.se", type = SE.type),SIMPLIFY = F)
        out <- lapply(out,function(x){rownames(x) <- c("Est.","S.E.");round(x,digits)})
     }else{
      out <- lapply(internalextract(object,what = "itemprob.parm"),round,digits)
    }
  }else if(tolower(what)=="lcprob"){

      out <- round(internalextract(object,what = "LC.prob"),digits)

}else if(tolower(what)=="gs"){
    if(withSE){
      gs <- t(sapply(internalextract(object,what = "catprob.parm"),
                      function(x)c(x[1],1-rev(x)[1])))
      se <- t(sapply(internalextract(object,what = "catprob.se", type = SE.type),
                     function(x)c(x[c(1,length(x))])))
      out <- round(cbind(gs,se),digits)
      colnames(out) <- c("guessing","slip","SE[guessing]","SE[slip]")
    }else{
    out <- round(t(sapply(internalextract(object,what = "catprob.parm"),
                    function(x)c(x[1],1-rev(x)[1]))),digits)
    colnames(out) <- c("guessing","slip")
    }
  }else if(tolower(what)=="delta"){
    d <- format_delta(delta = internalextract(object,what = "delta.parm"),
                      model = internalextract(object,what = "models_numeric"),
                      Kj = internalextract(object,what = "Kj"),
                      item.names = internalextract(object,what = "item.names"),
                      digits = digits+1)
    if(withSE){
      out <- mapply(rbind,d,
                    internalextract(object,what = "delta.se", type = SE.type),SIMPLIFY = F)
      out <- lapply(out,function(x){rownames(x) <- c("Est.","S.E.");round(x,digits)})
    }else{
      out <- lapply(d,round,digits)
    }
  }else if(tolower(what)=="rrum"){
    J <- internalextract(object,what = "ncat")
    ip <- internalextract(object,what = "catprob.parm")
    models <- internalextract(object,what = "models_numeric")
    if(all(models!=5)) stop("RRUM parameters are only available when RRUM is fitted.",call. = FALSE)
    out <- vector("list",J)
      for (j in 1:J){
        if (models[j]==5){
          parj <- ip[[j]]
          Lj <- length(parj)
          Kj <- log(Lj,2)
          pi.star <- parj[Lj]
          out[[j]] <- round(c(pi.star,rev(parj[(Lj-Kj):(Lj-1)])/pi.star),digits)
          names(out[[j]]) <- c("pi*",paste0("r",1:Kj))
        }

      }
      names(out) <- internalextract(object,what = "item.names")
      if(withSE) message("Standard errors are not available for RRUM parameters.")
    }else{
      stop(sprintf("No item parameters called \'%s\'", what), call.=FALSE)
  }

  return(out)

}
