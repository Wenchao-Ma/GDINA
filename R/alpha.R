#' Generate all possible attribute patterns
#'
#' This function generates all possible attribute patterns.
#' The Q-matrix is needed when any attributes are polytomous.
#'
#' @param K number of attributes
#' @param poly logical; is Q matrix polytomous?
#' @param Q Q-matrix; required when Q-matrix is polytomous
#'
#' @return attribute profiles for \eqn{2^K} latent classes
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @export
#' @examples
#' attributepattern(3)
#'
#' q <- matrix(scan(text = "0 1 2 1 0 1 1 2 0"),ncol = 3)
#' q
#' attributepattern(ncol(q),poly=TRUE,q)
#'
#' q <- matrix(scan(text = "0 1 1 1 0 1 1 1 0"),ncol = 3)
#' q
#' attributepattern(ncol(q),poly=TRUE,q)

attributepattern <-  function(K,poly=F,Q=NULL){
  return(alpha(K,poly,Q))
}



