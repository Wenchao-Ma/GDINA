#' Classification Rate Evaluation
#'
#' This function evaluates the classification rates for two sets of attribute profiles
#'
#' @param att1 a matrix or data frame of attribute profiles
#' @param att2 a matrix or data frame of attribute profiles
#'
#' @return a list with the following components:
#'  \describe{
  #' \item{PCA}{the proportion of correctly classified attributes (i.e., attribute level classification rate)}
  #' \item{PCV}{a vector giving the proportions of correctly classified attribute vectors (i.e., vector level classification rate).
  #'                   The fist element is the proportion of at least one attribute in the vector are correctly identified; the second
  #'                   element is the proportion of at least two attributes in the vector are correctly identified; and so forth. The last
  #'                   element is the proportion of all elements in the vector are correctly identified.}
  #' }
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu}}
#' @references
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive Diagnosis Modeling. \emph{Journal of Statistical Software, 93(14)}, 1-26.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' N <- 2000
#' # model does not matter if item parameter is probability of success
#' Q <- sim30GDINA$simQ
#' J <- nrow(Q)
#' gs <- matrix(0.1,J,2)
#'
#' set.seed(12345)
#' sim <- simGDINA(N,Q,gs.parm = gs)
#' GDINA.est <- GDINA(sim$dat,Q)
#'
#' CR <- ClassRate(sim$attribute,personparm(GDINA.est))
#' CR
#' }
ClassRate <- function(att1,att2){
  if (any(dim(att1)!=dim(att2))) stop("att1 and att2 must have the same dimensions.",call. = FALSE)
  comp <- as.matrix(att1)==as.matrix(att2)
  K <- ncol(att1)
  PCV <- NULL
  for (k in 1:K){
    PCV <- c(PCV,mean(rowSums(comp)>=k))
  }
  PCA <- mean(comp)
  return(list(PCA=PCA,PCV=PCV))
}
