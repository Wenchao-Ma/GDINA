#'@include GDINA-package.R GDINA.R
#'@title extract elements from objects of various classes
#'
#' @description A generic function to extract elements from objects of class \code{GDINA},
#' \code{itemfit}, \code{modelcomp}, \code{Qval} or \code{simGDINA}. This
#' page gives the elements that can be extracted from the class \code{GDINA}.
#' To see what can be extracted from \code{\link{itemfit}}, \code{\link{modelcomp}}, and
#' \code{\link{Qval}}, go to the corresponding function help page.
#'
#' Objects which can be extracted from \code{GDINA} objects include:
#'
#' \describe{
#' \item{itemprob.parm}{item success probability for each latent group}
#' \item{delta.parm}{delta parameters for each category}
#' \item{catprob.parm}{category success probability for each latent group; the same as itemprob.parm for dichotomous response items.}
#' \item{LCprob.parm}{category success probability for each latent class}
#' \item{itemprob.se}{SE associated with the item success probability for each latent group; need to specify either type=1 or type=2.}
#' \item{catprob.se}{SE associated with the category success probability for each latent group; need to specify either type=1 or type=2.}
#' \item{delta.se}{SE associated with the category success probability for each latent group; need to specify either type=1 or type=2.}
#' \item{higher.order.struc.parm}{higher-order structural parameters}
#' \item{higher.order.struc.se}{SE associated with higher-order structural parameters}
#'   \item{logLik}{observed log-likelihood value}
#'   \item{deviance}{deviance: -2 times observed log-likelihood value}
#'   \item{AIC}{AIC}
#'   \item{BIC}{BIC}
#'   \item{models}{fitted CDMs for each item/category}
#'   \item{npar}{number of parameters}
#'   \item{npar.item}{number of item parameters}
#'   \item{npar.att}{number of attribute parameters}
#'   \item{nobs}{number of individuals}
#'   \item{nitem}{number of items}
#'   \item{ncat}{number of categories excluding category zero}
#'   \item{natt}{number of attributes}
#'   \item{nitr}{number of iterations}
#'   \item{discrim}{GDINA discrimination index}
#'   \item{time}{time used}
#'   \item{start.time}{starting time}
#'   \item{end.time}{end time}
#'   \item{dat}{item responses analyzed}
#'   \item{Q}{Q-matrix}
#'   \item{Qc}{Qc-matrix}
#'   \item{posterior.prob}{posterior weights for each latent class}
#'   \item{prevalence}{prevalence of each attribute}
#' \item{itemprob.cov}{variance-covariance matrix of item endorsement probabilities for all items}
#' \item{itemprob.covindex}{index of the variance-covariance matrix of item endorsement probabilities for all items}
#' \item{logposterior.i}{log-posteriori for each examinee}
#' \item{loglikelihood.i}{log-likelihood for each examinee}
#' \item{expectedCorrect}{expected # of examinees in each latent group answering item correctly}
#' \item{expectedTotal}{expected # of examinees in each latent group}
#' \item{dif.LL}{absolute change in deviance in the last EM iteration}
#' \item{dif.p}{max absolute change in success probabilities in the last EM iteration}
#' \item{higher.order}{argument higher.order}
#' \item{higher.order.model}{argument higher.order.model}
#' \item{mono.constraint}{argument mono.constraint}
#' \item{SE}{argument SE}
#' \item{SE.type}{argument SE.type}
#' \item{empirical}{argument empirical}
#' \item{att.prior}{argument att.prior}
#' \item{att.str}{argument att.str}
#' \item{nstarts}{argument nstarts}
#' \item{conv.crit}{argument conv.crit}
#' \item{maxitr}{argument maxitr}
#' \item{higher.order.method}{argument higher.order.method}
#' \item{verbose}{argument verbose}
#' \item{sequential}{argument sequential}
#' \item{higher.order.parm}{argument higher.order.parm}
#' \item{digits}{argument digits}
#' \item{item.names}{argument item.names}
#' \item{itemprob.history}{itemprob.history in diagnosis mode}
#' \item{RN.history}{RN.history in diagnosis mode}
#' \item{likepost.history}{likepost.history in diagnosis mode}
#' \item{iter.history}{iter.history in diagnosis mode}
#' \item{HO.parm.history}{HO.parm.history in diagnosis mode}
#' \item{call}{function call}
#' }
#' @param object objects from class \code{GDINA},\code{itemfit}, \code{modelcomp}, \code{Qval} or \code{simGDINA}
#' @param what what to extract
#' @param digits how many decimal places for the output?
#' @param ... additional arguments
#'
#' @export
extract <- function (object, what, digits = 4,...) {
  UseMethod("extract")
}
#' @title extract elements of GDINA estimates
#' @description
#' NULL
#' @details
#' NULL
#' @describeIn GDINA extract various elements of GDINA estimates
#' @aliases extract.GDINA
#' @export
extract.GDINA <- function(object,what,digits=4,...){
  out <- switch(what,
                itemprob.parm = lapply(internalextract(object,"itemprob.parm"),round,digits),
                catprob.parm = lapply(internalextract(object,"catprob.parm"),round,digits),
                LCprob.parm = round(object$LC.prob,digits),
                higher.order.struc.parm = round(internalextract(object,"higher.order.struc.parm"),digits),
                catprob.se = {
                  se <- lapply(internalextract(object,"catprob.se",...),round,digits)
                  names(se) <- object$options$item.names
                  se
                  },
                itemprob.se = {
                  se <- lapply(internalextract(object,"itemprob.se",...),round,digits)
                  names(se) <- object$options$item.names
                  se
                },
                delta.parm = format_delta(internalextract(object,"delta.parm"),model = object$options$model,
                                          Kj = rowSums(object$options$Q),
                                          item.names = object$options$item.names,
                                          digits = digits),
                delta.se = {
                  format_delta(delta = delta_se(object,...)$se,
                               model = object$options$model, Kj = rowSums(internalextract(object,"Q")),
                               item.names = object$options$item.names,digits = digits)
                  },
                logLik = round(-0.5*object$testfit$Deviance,2),
                deviance = object$testfit$Deviance,
                npar = object$options$npar,
                npar.item = object$options$item.npar,
                npar.att = object$options$npar - object$options$item.npar,
                nitr = object$options$itr,
                nobs = nrow(object$options$dat),
                nitem = ncol(object$options$dat),
                ncat = nrow(object$options$Q),
                natt = ifelse(object$options$sequential,ncol(object$options$Q)-2,ncol(object$options$Q)),
                AIC = object$testfit$AIC,
                BIC = object$testfit$BIC,
                models = object$model,
                discrim = {
                  wp <- t(internalextract(object,"LCprob.parm")) * object$posterior.prob  #L x J w*p matrix
                  gdi <- colSums(t((object$LC.prob - colSums(wp))^2) * object$posterior.prob) # vector of length J
                  dj <- sapply(object$itemprob.parm, function(x) x[length(x)]-x[1])
                  Discrim <- round(data.frame(dj = dj, GDI = gdi),digits)
                  rownames(Discrim) <- object$options$item.names
                  colnames(Discrim) <- c("P(1)-P(0)", "GDI")
                  Discrim
                },
                time = object$options$timeused,
                start.time = object$options$start.time,
                end.time = object$options$end.time,
                dat = object$options$dat,
                Qc = internalextract(object,what = "Qc"),
                Q = internalextract(object,what = "Q"),
                posterior.prob = round(object$posterior.prob,digits),
                prevalence = round(internalextract(object,what = "prevalence"),digits),
                catprob.cov = internalextract(object,what = "catprob.cov",...),
                delta.cov = internalextract(object,what = "delta.cov",...),
                logposterior.i = object$technicals$logposterior.i,
                loglikelihood.i = object$technicals$loglikelihood.i,
                expectedCorrect = object$technicals$expectedCorrect,
                expectedTotal = object$technicals$expectedTotal,
                higher.order = object$options$higher.order,
                higher.order.model = object$options$higher.order.model,
                mono.constraint = object$options$mono.constraint ,
                Mstep = object$options$Mstep,
                person.est = object$options$person.est,
                SE = object$options$SE,
                SE.type = object$options$SE.type,
                empirical = object$options$empirical,
                att.prior = object$options$att.prior,
                att.str=object$options$att.str,
                nstarts = object$options$nstarts,
                conv.crit = object$options$conv.crit,
                maxitr = object$options$maxitr,
                higher.order.method = object$options$higher.order.method,
                higher.order.SE = object$options$higher.order.SE,
                verbose = object$options$verbose,
                sequential = object$options$sequential,
                higher.order.parm = object$options$higher.order.parm,
                digits = object$options$digits,
                dif.LL = object$options$dif.LL,
                dif.p=object$options$dif.p,
                item.names = object$options$item.names,
                itemprob.history = object$diagnos$itemprob.matrix,
                RN.history = object$diagnos$RN,
                likepost.history=object$diagnos$likepost,
                iter.history = object$diagnos$changelog,
                HO.parm.history = object$diagnos$HO.parm,
                call = object$options$call,
                expectedCorrect.LC = {
                  if(object$options$sequential){
                    out <- apply(seq_coding(object$options$dat,object$options$Q),2,
                          function(x)colSums(x*exp(object$technicals$logposterior.i),na.rm = TRUE))
                  }else{
                    out <- apply(object$options$dat,2,
                                 function(x)colSums(x*exp(object$technicals$logposterior.i),na.rm = TRUE))
                  }
                  out <- t(out)
                  row.names(out) <- object$options$item.names
                  out
                },
                expectedTotal.LC = {
                  if(object$options$sequential){
                    out <- apply(seq_coding(object$options$dat,object$options$Q),2,
                                 function(x)colSums(as.numeric(!is.na(x))*exp(object$technicals$logposterior.i),na.rm = TRUE))
                  }else{
                    out <- matrix(colSums(exp(object$technicals$logposterior.i),na.rm = TRUE),
                                  ncol = nrow(object$options$Q), nrow = ncol(object$technicals$logposterior.i))
                  }
                  out <- t(out)
                  row.names(out) <- object$options$item.names
                  colnames(out) <- colnames(object$technicals$logposterior.i)
                  out
                },
                stop(sprintf("Can not extract element \'%s\'", what), call.=FALSE))
  return(out)
}
