#' extract various components from model comparison
#'
#' @description
#'NULL
#'
#'
#' @param object object of class \code{modelcomp} for various S3 methods
#' @param what argument for S3 method \code{extract} indicating what to extract;
#' It can be \code{"wald"} for wald statistics, \code{"wald.p"} for associated p-values,
#' \code{"df"} for degrees of freedom,
#' and \code{"DS"} for dissimilarity between G-DINA and other CDMs.
#' @param digits How many decimal places in each number? The default is 4.
#' @param ... additional arguments
#'
#' @author Wenchao Ma & Jimmy de la Torre
#' @include modelcomp.R
#'@describeIn modelcomp extract various elements from \code{modelcomp} objects
#'@aliases extract.modelcomp
#'@export
extract.modelcomp <- function(object,
                           what=c("wald","wald.p","df","DS","models"), digits = 4,...){
  what <- match.arg(what)
  out <- switch(what,
                wald = round(object$wald,digits),
                wald.p = round(object$wald.p,digits),
                df = object$df,
                DS = {if(is.null(object$DS)) NULL else round(object$DS,digits)},
                models = object$models)
  return(out)
}
