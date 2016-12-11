#'@title NULL
#'
#'@description NULL
#'
#' @param object \code{Qval} objects for S3 methods
#' @param what argument for S3 method \code{extract} indicating what to extract;
#' It can be \code{"sug.Q"} for suggested Q-matrix,
#' \code{"Q"} for original Q-matrix, \code{"varsigma"} for varsigma index,
#' and \code{"PVAF"} for PVAF.
#' @param ... additional arguments
#'
#'@describeIn Qval extract various elements from \code{Qval} objects
#'@aliases extract.Qval
#'@export
extract.Qval <- function(object,
                      what=c("sug.Q","varsigma","PVAF","eps","Q"),...){
  what <- match.arg(what)
  out <- switch(what,
                sug.Q = object$sug.Q,
                varsigma = object$varsigma,
                PVAF = object$PVAF,
                eps = object$eps,
                Q = object$Q)
  return(out)
}
