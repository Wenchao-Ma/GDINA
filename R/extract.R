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
#' @include modelcomp.R
#'@describeIn modelcomp extract various elements from \code{modelcomp} objects
#'@aliases extract.modelcomp
#'@export
extract.modelcomp <- function(object,
                              what=c("stats","pvalues","adj.pvalues","df","DS","models"), digits = 4,...){
  what <- match.arg(what)
  out <- switch(what,
                stats = round(object$stats,digits),
                pvalues = round(object$pvalues,digits),
                adj.pvalues = round(object$adj.pvalues,digits),
                df = object$df,
                DS = {if(is.null(object$DS)) NULL else round(object$DS,digits)},
                models = object$models)
  return(out)
}



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
                Q = object$Q,
                stop(sprintf("Can not show element \'%s\'", what), call.=FALSE))
  return(out)
}




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
#'   \item{att.prior}{attribute prior weights for calculating marginalized likelihood in the last EM iteration}
#'   \item{discrim}{GDINA discrimination index}
#'   \item{designmatrix}{A list of design matrices for each item/category}
#' \item{expectedCorrect}{expected # of examinees in each latent group answering item correctly}
#' \item{expectedTotal}{expected # of examinees in each latent group}
#' \item{higher.order}{higher-order model specifications}
#' \item{linkfunc}{link functions for each item}
#' \item{initial.catprob}{initial item category probability parameters}
#'   \item{prevalence}{prevalence of each attribute}
#'   \item{posterior.prob}{posterior weights for each latent class}
#' }
#' @param object objects from class \code{GDINA},\code{itemfit}, \code{modelcomp}, \code{Qval} or \code{simGDINA}
#' @param what what to extract
#' @param ... additional arguments
#'
#' @examples
#'\dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' extract(fit,"discrim")
#' extract(fit,"designmatrix")
#' }
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
                  if(SE.type==3&any(extract(object,"att.dist")=="higher.order"))stop("standard error cannot be calculated.",call. = FALSE)
                  Kj <- extract(object, "Kj")
                  se <- OPG_p(object, SE.type = SE.type, ...)$se
                  for (j in 1:length(se)) names(se[[j]]) <- paste0("SE[P(", apply(attributepattern(Kj[j]), 1, paste0, collapse = ""), ")]")
                  names(se) <- object$options$item.names
                  se
                },
                dat = object$options$dat, #raw data
                delta.cov = {var <- OPG_d(object, SE.type = SE.type);list(cov = var$cov, index = var$ind)},
                delta.parm = {format_delta(object$delta.parm,model = object$options$model,
                                           Kj = rowSums(extract(object,"Q")),
                                           item.names = object$options$item.names,
                                           digits = 10)},
                delta.se = {
                  if(SE.type==3&any(extract(object,"att.dist")=="higher.order"))stop("standard error cannot be calculated.",call. = FALSE)
                  format_delta(delta = OPG_d(object,SE.type = SE.type)$se,
                               model = object$options$model, Kj = rowSums(extract(object,"Q")),
                               item.names = object$options$item.names,digits = 10)
                },
                designmatrix = object$options$DesignMatrices,
                deviance = object$testfit$Deviance,
                dif.p = object$options$dif.p,
                dif.prior = object$options$dif.prior,
                dif.LL = object$options$dif.LL,
                discrim = {
                  if(object$options$no.group>1) stop("Discrimination indices are not available for multiple group model.",call. = FALSE)
                  dj <- sapply(object$catprob.parm, function(x) x[length(x)]-x[1])
                  wp <- t(object$LC.prob) * c(object$posterior.prob)  #L x J w*p matrix
                  gdi <- colSums(t((object$LC.prob - colSums(wp))^2) * c(object$posterior.prob)) # vector of length J
                  Discrim <- data.frame(dj = dj, GDI = gdi)
                  rownames(Discrim) <- object$options$item.names
                  colnames(Discrim) <- c("P(1)-P(0)", "GDI")
                  Discrim
                },
                end.time = object$extra$end.time,
                expectedCorrect = object$technicals$expectedCorrect,
                expectedTotal = object$technicals$expectedTotal,
                expectedCorrect.LC = {
                  if(object$options$sequential){
                    dat <- extract(object,"seq.dat")
                  }else{
                    dat <- extract(object,"dat")
                  }
                  out <- NgRg(as.matrix(object$technicals$logposterior.i),
                              as.matrix(dat),
                              eta(matrix(1,extract(object,"ncat"),extract(object,"natt"))),
                              rep(1,nrow(dat)))$Rg
                  row.names(out) <- object$options$item.names
                  out
                },
                expectedTotal.LC = {
                  if(object$options$sequential){
                    dat <- extract(object,"seq.dat")
                  }else{
                    dat <- extract(object,"dat")
                  }
                  out <- NgRg(as.matrix(object$technicals$logposterior.i),
                              as.matrix(dat),
                              eta(matrix(1,extract(object,"ncat"),extract(object,"natt"))),
                              rep(1,nrow(dat)))$Ng
                  row.names(out) <- object$options$item.names
                  out
                },
                group = object$options$group,
                gr = object$options$gr,
                higher.order = object$options$higher.order,
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
                          pj <- sj <- rbind(uP(as.matrix(eta(as.matrix(redQj))), as.matrix(itemparj)),0)
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
                      colnames(pj) <- paste0("P(", apply(attributepattern(Kj), 1, paste0, collapse = ""), ")")
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
                Kj = {rowSums(extract(object,"Q"))},
                latent.var = object$options$latent.var,
                LCprob.parm = object$LC.prob,
                LCpf.parm = {
                  LCpf <- patt <- eta(as.matrix(extract(object,"Q")))
                  pf <- extract(object,"catprob.parm")
                  for(i in 1:nrow(LCpf)) LCpf[i,] <- pf[[i]][patt[i,]]
                  rownames(LCpf) <- rownames(extract(object,"catprob.matrix"))
                  LCpf
                },# processing function
                linkfunc = {c("identity","logit","log")[object$options$linkfunc]},
                logLik = -0.5*object$testfit$Deviance,
                logposterior.i = object$technicals$logposterior.i,
                loglikelihood.i = object$technicals$loglikelihood.i,
                loglinear = object$options$loglinear,
                models = object$model,
                models_numeric = object$options$model,
                mono.constraint = object$options$mono.constraint,
                npar = object$technicals$total.npar,
                npar.item = object$technicals$free.item.npar, # free parameters
                npar.fixeditem = object$technicals$total.item.npar - object$technicals$free.item.npar,
                npar.att = object$technicals$stru.npar,
                natt = ncol(extract(object,"Q")),
                ncat = nrow(extract(object,"Q")),
                ngroup = object$options$no.group,
                nitem = ncol(object$options$dat),
                nitr = object$options$itr,
                nobs = nrow(object$options$dat),
                pf = object$pf,
                posterior.prob = object$posterior.prob,
                prevalence = {
                  Q <- extract(object,"Q")
                  preva <- vector("list",extract(object,"ngroup"))
                  pattern <- attributepattern(Q = Q)
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
                  if (object$options$sequential||any(extract(object,"models")=="MSDINA")) {
                    out <- data.frame(object$options$Qm)
                  } else{
                    out <- data.frame(object$options$Q)
                  }
                  colnames(out) <- paste0("A", 1:ncol(out))
                  out
                },
                Qc = {if(any(extract(object,"models")=="MSDINA")){
                  object$options$Qcm
                } else if(object$options$sequential) {
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
                originalQ = object$options$Q,
                start.time = object$extra$start.time,
                time = object$extra$timeused,
                sequential = object$options$sequential,
                seq.dat = object$options$seq.dat,
                struc.parm = object$struc.parm,
                verbose = object$options$verbose,
                call = object$extra$call,
                stop(sprintf("Can not extract element \'%s\'", what), call.=FALSE))
  return(out)
}
