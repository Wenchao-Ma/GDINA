#' @include simGDINA.R
#' @title NULL
#'
#' @description  NULL
#'
#' @param object object of class \code{simGDINA} for method \code{extract}
#' @param what argument for S3 method \code{extract} indicating what to extract
#' @param ... additional arguments
#' @rdname simGDINA
#'@export
extract.simGDINA <- function(object,
                          what=c("dat","Q","attribute","catprob.parm",
                                 "delta.parm","higher.order.parm","mvnorm.parm",
                                 "LCprob.parm"),...){
  out <- switch(tolower(what),
                dat = object$dat,
                Q = object$Q,
                attribute = object$attribute,
                catprob.parm = object$catprob.parm,
                delta.parm = object$delta.parm,
                higher.order.parm = object$higher.order.parm,
                mvnorm.parm = object$mvnorm.parm,
                LCprob.parm = object$LCprob.parm,
                stop(sprintf("Can not show element \'%s\'", what), call.=FALSE))
  return(out)
}
