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
#'   \item{AIC}{AIC}
#'   \item{att.prior}{attribute prior weights for calculating marginalized likelihood in the last iteration}
#'  \item{att.str}{argument att.str}
#'   \item{BIC}{BIC}
#' \item{call}{function call}
#' \item{catprob.parm}{category success probability for each latent group; the same as itemprob.parm for dichotomous response items.}
#' \item{catprob.se}{SE associated with the category success probability for each latent group.}
#' \item{catprob.cov}{variance-covariance matrix of item endorsement probabilities for all items}
#' \item{conv.crit}{argument conv.crit}
#'  \item{dat}{item responses analyzed}
#' \item{delta.parm}{delta parameters for each category}
#' \item{delta.cov}{Convariance matrix associated with the delta parameters.}
#' \item{delta.se}{SE associated with the delta parameters for each latent group.}
#'   \item{deviance}{deviance: -2 times observed log-likelihood value}
#' \item{dif.LL}{absolute change in deviance in the last EM iteration}
#' \item{dif.p}{max absolute change in success probabilities in the last EM iteration}
#' \item{digits}{argument digits}
#'   \item{discrim}{GDINA discrimination index}
#' \item{empirical}{argument empirical}
#'   \item{end.time}{end time}
#' \item{expectedCorrect}{expected # of examinees in each latent group answering item correctly}
#' \item{expectedTotal}{expected # of examinees in each latent group}
#' \item{higher.order}{higher-order model specifications}
#' \item{higher.order.method}{argument higher.order$method}
#' \item{higher.order.model}{argument higher.order$model}
#' \item{HO.parm.history}{HO.parm.history in diagnosis mode}
#' \item{initial.catprob}{initial item category probability parameters}
#' \item{iter.history}{iter.history in diagnosis mode}
#' \item{item.names}{argument item.names}
#' \item{itemprob.history}{itemprob.history in diagnosis mode}
#' \item{itemprob.parm}{item success probability for each latent group}
#' \item{itemprob.se}{SE associated with the item success probability for each latent group.}
#' \item{LCprob.parm}{category success probability for each latent class}
#'   \item{logLik}{observed log-likelihood value}
#' \item{loglikelihood.i}{log-likelihood for each examinee}
#' \item{likepost.history}{likepost.history in diagnosis mode}
#' \item{logposterior.i}{log-posteriori for each examinee}
#' \item{maxitr}{argument maxitr}
#'   \item{models}{fitted CDMs for each item/category}
#' \item{mono.constraint}{argument mono.constraint}
#'   \item{natt}{number of attributes}
#'   \item{ncat}{number of categories excluding category zero}
#'   \item{ngroup}{number of groups}
#'   \item{nitem}{number of items}
#'   \item{nitr}{number of iterations}
#'   \item{nobs}{number of individuals}
#'   \item{npar}{number of parameters}
#'   \item{npar.item}{number of item parameters}
#'   \item{npar.att}{number of attribute parameters}
#' \item{nstarts}{argument nstarts}
#'   \item{prevalence}{prevalence of each attribute}
#'   \item{posterior.prob}{posterior weights for each latent class}
#'   \item{Q}{Q-matrix}
#'   \item{Qc}{Qc-matrix}
#' \item{RN.history}{RN.history in diagnosis mode}
#'   \item{start.time}{starting time}
#' \item{sequential}{argument sequential}
#' \item{seq.dat}{data for sequential models}
#'   \item{time}{time used}
#' \item{verbose}{argument verbose}
#' }
#' @param object objects from class \code{GDINA},\code{itemfit}, \code{modelcomp}, \code{Qval} or \code{simGDINA}
#' @param what what to extract
#' @param ... additional arguments
#'
#' @export
extract <- function (object, what, ...) {
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
extract.GDINA <- function(object,what,SE.type = 2,...){
  out <- switch(what,
                AIC = object$testfit$AIC,
                att.dist = object$options$att.dist,
                att.prior = object$options$att.prior,
                att.str=object$options$att.str,
                BIC = object$testfit$BIC,
                catprob.cov = {var <- OPG_p(object, SE.type = SE.type);list(cov = var$cov, index = var$ind)},
                catprob.matrix = object$catprob.matrix,
                catprob.parm = object$catprob.parm,
                catprob.se = {
                  Kj <- extract(object, "Kj")
                  se <- OPG_p(object, SE.type = SE.type, ...)$se
                  for (j in 1:length(se)) names(se[[j]]) <- paste0("SE[P(", apply(alpha(Kj[j]), 1, paste0, collapse = ""), ")]")
                  names(se) <- object$options$item.names
                  se
                },
                conv.crit = object$options$conv.crit,
                dat = object$options$dat,
                delta.cov = {var <- OPG_d(object, SE.type = SE.type);list(cov = var$cov, index = var$ind)},
                delta.parm = {format_delta(object$delta.parm,model = object$options$model,
                                          Kj = rowSums(extract(object,"Q")),
                                          item.names = object$options$item.names,
                                          digits = 10)},
                delta.se = {
                  format_delta(delta = OPG_d(object,SE.type = SE.type)$se,
                               model = object$options$model, Kj = rowSums(extract(object,"Q")),
                               item.names = object$options$item.names,digits = 10)
                },
                deviance = object$testfit$Deviance,
                digits = object$options$digits,
                dif.p=object$options$dif.p,
                dif.LL = object$options$dif.LL,
                discrim = {
                  if(object$options$no.group>1) warning("Discrimination indices are not available for multiple group model.",call. = FALSE);return(NULL)
                  dj <- sapply(object$catprob.parm, function(x) x[length(x)]-x[1])
                  wp <- t(object$LC.prob) * c(object$posterior.prob)  #L x J w*p matrix
                  gdi <- colSums(t((object$LC.prob - colSums(wp))^2) * c(object$posterior.prob)) # vector of length J
                  Discrim <- data.frame(dj = dj, GDI = gdi)
                  rownames(Discrim) <- object$options$item.names
                  colnames(Discrim) <- c("P(1)-P(0)", "GDI")
                  Discrim
                },
                empirical = object$options$empirical,
                end.time = object$options$end.time,
                expectedCorrect = object$technicals$expectedCorrect,
                expectedTotal = object$technicals$expectedTotal,
                expectedCorrect.LC = {
                if(object$options$sequential){
                    dat <- extract(object,"seq.dat")
                    }else{
                      dat <- extract(object,"dat")
                    }
                  out <- NgRg(object$technicals$logposterior.i,dat,
                       eta.loc(matrix(1,extract(object,"ncat"),extract(object,"natt"))))$Rg
                  row.names(out) <- object$options$item.names
                  out
                },
                expectedTotal.LC = {
                  # if(object$options$sequential){
                  #   out <- apply(seq_coding(object$options$dat,object$options$Q),2,
                  #                function(x)colSums(as.numeric(!is.na(x))*exp(object$technicals$logposterior.i),na.rm = TRUE))
                  # }else{
                  #   out <- matrix(colSums(exp(object$technicals$logposterior.i),na.rm = TRUE),
                  #                 ncol = nrow(object$options$Q), nrow = ncol(object$technicals$logposterior.i))
                  # }
                  # out <- t(out)
                  if(object$options$sequential){
                    dat <- extract(object,"seq.dat")
                  }else{
                    dat <- extract(object,"dat")
                  }
                  out <- NgRg(object$technicals$logposterior.i,dat,
                              eta.loc(matrix(1,extract(object,"ncat"),extract(object,"natt"))))$Ng
                  row.names(out) <- object$options$item.names
                  # colnames(out) <- colnames(object$technicals$logposterior.i)
                  out
                },
                higher.order = object$options$higher.order,
                higher.order.model = object$options$higher.order.model,
                higher.order.method = object$options$higher.order.method,
                higher.order.SE = object$options$higher.order.SE,
                higher.order.struc.parm = object$higher.order.struc.parm,
                initial.catprob = object$technical$initial.parm,
                itemprob.parm = {
                  if(object$options$sequential){
                    itemparm <- object$catprob.matrix
                    p <- vector("list",extract(object,"nitem"))
                    Qc <- extract(object,"Qc")
                    for (j in 1:extract(object,"nitem")){
                      Qj <- Qc[which(Qc[,1]==j),-c(1:2),drop=FALSE]
                      itemparj <- itemparm[which(Q[,1]==j),,drop=FALSE]
                      colj <- which(colSums(Qj)>0)
                      Kj <- length(colj)
                      if(nrow(Qj)>1){ # polytomous items
                        if(length(colj)>1){
                          redQj <- Qj[,colj]
                          pj <- sj <- rbind(uP(as.matrix(eta.loc(redQj)), as.matrix(itemparj)),0)
                        }else{
                          pj <- sj <- rbind(itemparj[,1:2],0)
                        }
                        for (s in 1:(nrow(sj)-1)){
                          pj[s,] <- apply(sj[1:s,,drop=FALSE],2,prod)*(1-sj[s+1,])
                        }
                        pj <- pj[-nrow(pj),]
                      }else{ #dichotomous items
                        pj <- matrix(itemparj[!is.na(itemparj)],nrow = 1)
                      }
                      colnames(pj) <- paste0("P(", apply(alpha(Kj), 1, paste0, collapse = ""), ")")
                      rownames(pj) <- paste("Cat",1:nrow(pj))
                      p[[j]] <- pj
                    }
                    names(p) <- paste("Item",1:extract(object,"nitem"))
                    p
                  }else{
                    object$catprob.parm
                  }

                },
                itemprob.se = { if(!object$options$sequential) extract(object,"catprob.se",SE.type = SE.type)},
                item.names = object$options$item.names,
                itemprob.history = object$diagnos$itemprob.matrix,
                iter.history = object$diagnos$changelog,
                Kj = {rowSums(extract(object,"Q"))},
                LCprob.parm = object$LC.prob,
                logLik = -0.5*object$testfit$Deviance,
                logposterior.i = object$technicals$logposterior.i,
                loglikelihood.i = object$technicals$loglikelihood.i,
                likepost.history=object$diagnos$likepost,
                maxitr = object$options$maxitr,
                models = object$model,
                models_numeric = object$options$model,
                mono.constraint = object$options$mono.constraint ,
                Mstep = object$options$Mstep,
                npar = object$options$npar,
                npar.item = object$options$item.npar,
                npar.att = object$options$npar - object$options$item.npar,
                natt = ncol(extract(object,"Q")),
                ncat = nrow(extract(object,"Q")),
                ngroup = object$options$no.group,
                nitem = ncol(object$options$dat),
                nitr = object$options$itr,
                nobs = nrow(object$options$dat),
                nstarts = object$options$nstarts,
                posterior.prob = object$posterior.prob,
                prevalence = {
                  Q <- extract(object,"Q")
                  preva <- vector("list",extract(object,"ngroup"))
                  pattern <- alpha(ncol(Q), T, Q)
                  for(g in 1:extract(object,"ngroup")) {
                    for (i in c(0:max(Q))) {
                      preva[[g]] <-
                        cbind(preva[[g]], t(extract(object,"posterior.prob")[g,,drop=FALSE] %*% ((pattern == i) * 1)))
                    }
                    colnames(preva[[g]]) <- paste0("Level", 0:max(Q))
                    rownames(preva[[g]]) <- paste0("A", seq(1, ncol(Q)))
                  }
                  names(preva) <- object$options$group.label
                  preva
                },
                Q = {
                  if (object$options$sequential) {
                    out <- data.frame(object$options$Q[, -c(1:2)])
                  } else{
                    out <- data.frame(object$options$Q)
                  }
                  colnames(out) <- paste0("A", 1:ncol(out))
                  out
                },
                Qc = {if (object$options$sequential) {
                  data.frame(object$options$Q)
                } else{
                  out <- data.frame(
                    Item = 1:nrow(object$options$Q),
                    Cat = rep(1, nrow(object$options$Q)),
                    object$options$Q
                  )
                  colnames(out)[-c(1:2)] <-
                    paste0("A", 1:ncol(object$options$Q))
                  out
                }},
                RN.history = object$diagnos$RN,
                start.time = object$options$start.time,
                time = object$options$timeused,
                sequential = object$options$sequential,
                seq.dat = object$options$seq.dat,
                verbose = object$options$verbose,
                HO.parm.history = object$diagnos$HO.parm,
                call = object$options$call,
                stop(sprintf("Can not extract element \'%s\'", what), call.=FALSE))
  return(out)
}
