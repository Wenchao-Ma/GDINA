#'@include itemfit.R
#' extract elements from objects of itemfit
#'
#'@description NULL
#'
#' @param object objects of class \code{itemfit} for various S3 methods
#' @param what argument for S3 method \code{extract} indicating what to extract;
#' It can be \code{"p"} for proportion correct statistics,
#' \code{"r"} for transformed correlations, \code{logOR} for log odds ratios and
#' \code{"maxitemfit"} for maximum statistics for each item.
#' @param ... additional arguments
#'
#' @describeIn itemfit extract various elements from \code{itemfit} objects
#' @aliases extract.itemfit
#' @export
extract.itemfit <- function(object,what,...){
  out <- switch(what,
                r = object$r,
                p = object$p,
                logOR = object$logOR,
                maxitemfit = object$max.itemlevel.fit,
                person.sim = object$options$person.sim,
                p.adjust.method = object$options$p.adjust.methods,
                N.resampling = object$options$N.resampling,
                randomseed = object$options$randomseed,
                call = object$options$call,
                digits = object$options$digits,
                stop(sprintf("Can not extract element \'%s\'", what), call.=FALSE))
  return(out)
}
