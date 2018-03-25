#'@title extract item parameters (deprecated)
#'
#' @description
#'
#' This function has been deprecated; use \code{coef} instead.
#'
#' @param object estimated GDINA object returned from \code{\link{GDINA}}
#' @param what what to show.
#' @param withSE show standard errors or not?
#' @param SE.type Type of standard errors.
#' @param digits how many decimal places for the ouput?
#' @param ... additional arguments
#'
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' # deprecated
#' itemparm(fit)
#' coef(fit)
#'}
#'
#' @references
#'
#' Philipp, M., Strobl, C., de la Torre, J., & Zeileis, A.(2017). On the estimation of standard errors in cognitive diagnosis models. \emph{Journal of Educational and Behavioral Statistics, 43}, 88-115.
#'
#'@export
itemparm <- function(object,
                      what = c("catprob","gs","delta","rrum","itemprob","LCprob"),
                      withSE = FALSE, SE.type = 2,digits = 4, ...) {
  UseMethod("itemparm")
}


#' @rdname itemparm
#' @export
itemparm.GDINA <- function(object,
                           what = c("catprob","gs","delta","rrum","itemprob","LCprob"),
                           withSE = FALSE, SE.type = 2,digits = 4, ...){
  if(!class(object)=="GDINA") stop("object must be a GDINA estimate.",call. = FALSE)
  .Deprecated("coef", package="GDINA",msg = "'itemparm' is deprecated - use 'coef' instead.")
  what <- match.arg(what)
  if(tolower(what)=="catprob"){
    if(withSE){
      out <- mapply(rbind,extract(object,what = "catprob.parm"),
                    extract(object,what = "catprob.se", SE.type = SE.type),SIMPLIFY = F)
      out <- lapply(out,function(x){rownames(x) <- c("Est.","S.E.");round(x,digits)})
    }else{
      out <- lapply(extract(object,what = "catprob.parm"),round,digits)
    }
  }else if(tolower(what)=="itemprob"){
    if(withSE&!extract(object,what = "sequential")){
        out <- mapply(rbind,extract(object,what = "itemprob.parm"),
                      extract(object,what = "itemprob.se", SE.type = SE.type),SIMPLIFY = F)
        out <- lapply(out,function(x){rownames(x) <- c("Est.","S.E.");round(x,digits)})
     }else{
      out <- lapply(extract(object,what = "itemprob.parm"),round,digits)
    }
  }else if(tolower(what)=="lcprob"){

      out <- round(extract(object,what = "LCprob.parm"),digits)

}else if(tolower(what)=="gs"){
    if(withSE){
      gs <- t(sapply(extract(object,what = "catprob.parm"),
                      function(x)c(x[1],1-rev(x)[1])))
      se <- t(sapply(extract(object,what = "catprob.se", SE.type = SE.type),
                     function(x)c(x[c(1,length(x))])))
      out <- round(cbind(gs,se),digits)
      colnames(out) <- c("guessing","slip","SE[guessing]","SE[slip]")
    }else{
    out <- round(t(sapply(extract(object,what = "catprob.parm"),
                    function(x)c(x[1],1-rev(x)[1]))),digits)
    colnames(out) <- c("guessing","slip")
    }
  }else if(tolower(what)=="delta"){
    d <- format_delta(delta = extract(object,what = "delta.parm"),
                      model = extract(object,what = "models_numeric"),
                      Kj = extract(object,what = "Kj"),
                      item.names = extract(object,what = "item.names"),
                      digits = digits+1)
    if(withSE){
      out <- mapply(rbind,d,
                    extract(object,what = "delta.se", SE.type = SE.type),SIMPLIFY = F)
      out <- lapply(out,function(x){rownames(x) <- c("Est.","S.E.");round(x,digits)})
    }else{
      out <- lapply(d,round,digits)
    }
  }else if(tolower(what)=="rrum"){
    J <- extract(object,what = "ncat")
    ip <- extract(object,what = "catprob.parm")
    models <- extract(object,what = "models_numeric")
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
      names(out) <- extract(object,what = "item.names")
      if(withSE) message("Standard errors are not available for RRUM parameters.")
    }else{
      stop(sprintf("No item parameters called \'%s\'", what), call.=FALSE)
  }

  return(out)

}
