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
#' @param digits how many decimal places for the output?
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
#' }
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
  stopifnot(isa(object,"GDINA"))
  .Deprecated("coef", package="GDINA",msg = "'itemparm' is deprecated - use 'coef' instead.")
  # Delegate to coef.GDINA after mapping deprecated parameter names
  coef(object, what = what, withSE = withSE, SE.type = SE.type, digits = digits, ...)
}
