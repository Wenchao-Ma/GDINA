#'@title Calculate classification accuracy
#'
#' @description
#' This function calculate test-, pattern- and attribute-level classification accuracy indices based on GDINA estimates from
#' the \code{GDINA} function using approaches in Iaconangelo (2017) and Wang, Song, Chen, Meng, and Ding (2015).
#' It is only applicable for dichotomous attributes.
#'
#' @param GDINA.obj estimated GDINA object returned from \code{\link{GDINA}}
#' @param what what attribute estimates are used? Default is \code{"MAP"}.
#'
#' @return a list with elements
#' \describe{
#' \item{tau}{estimated test-level classification accuracy, see Iaconangelo (2017, Eq 2.2)}
#' \item{tau_l}{estimated pattern-level classification accuracy, see Iaconangelo (2017, p. 13)}
#' \item{tau_k}{estimated attribute-level classification accuracy, see Wang, et al (2015, p. 461 Eq 6)}
#' \item{CCM}{Conditional classification matrix, see Iaconangelo (2017, p. 13)}
#' }
#'
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#' @references
#'
#' Iaconangelo, C.(2017). \emph{Uses of Classification Error Probabilities in the Three-Step Approach to Estimating Cognitive Diagnosis Models.} (Unpublished doctoral dissertation). New Brunswick, NJ: Rutgers University.
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive Diagnosis Modeling. \emph{Journal of Statistical Software, 93(14)}, 1-26.
#'
#' Wang, W., Song, L., Chen, P., Meng, Y., & Ding, S. (2015). Attribute-Level and Pattern-Level Classification Consistency and Accuracy Indices for Cognitive Diagnostic Assessment.
#' \emph{Journal of Educational Measurement, 52} , 457-476.
#'
#'
#' @examples
#'\dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' fit
#' CA(fit)
#' }
#'@export




CA <- function(GDINA.obj,what="MAP"){
  if(extract(GDINA.obj, "ngroup")!=1) {
    stop("classification accuracy is only applicable for single group analysis.", call. = FALSE)
  }
  p_c <- extract(GDINA.obj,"posterior.prob")
  pp <- personparm(GDINA.obj,what = what)
  if(what=="MAP"||what=="MLE"){
    if(any(pp[,ncol(pp)])) warning(paste0(what," estimates for some individuals have multiple modes.",collapse = ""),call. = FALSE)
    pp <- as.matrix(pp[,-ncol(pp)])
  }
  mp <- personparm(GDINA.obj,what = "mp")
  patt <- extract(GDINA.obj,"attributepattern")
  gr <- matchMatrix(patt,pp)
  pseudo.gr <- setdiff(seq(nrow(patt)),unique(gr))
  gr <- c(gr,pseudo.gr)
  lab <- apply(patt,1,paste0,collapse="")
  # conditional classification matrix
  # row: true; col: estimated
  post <- cbind(exp(t(indlogPost(GDINA.obj))),matrix(0,nrow(patt),length(pseudo.gr)))
  CCM <- aggregateCol(post,gr)/c(extract(GDINA.obj,"nobs")*p_c)
  tau_c <- diag(CCM)
  tau <- sum(tau_c*c(p_c))
  tau_k <- colMeans(pp*mp+(1-pp)*(1-mp))
  names(tau_c) <- rownames(CCM) <- colnames(CCM) <- lab
  ret <- list(tau=tau,tau_l=tau_c,tau_k=tau_k,CCM = CCM)
  class(ret) <- "CA"
  return(ret)
}
