#' @include GDINA.R
#' @title  Differential item functioning for cognitive diagnosis models
#'
#' @description   This function is used to detect differential item functioning based on the models estimated
#' in the \code{\link{GDINA}} function using the Wald test (Hou, de la Torre, & Nandakumar, 2014) and the likelihood ratio
#' test (Ma, Terzi, Lee, & de la Torre, 2017). It can only detect DIF for two groups currently.
#'
#' @param dat item responses from two groups; missing data need to be coded as \code{NA}
#' @param Q Q-matrix specifying the association between items and attributes
#' @param group a numerical vector with integer 1, 2, ..., # of groups indicating the group each individual belongs to. It must start from 1 and its
#'    length must be equal to the number of individuals.
#' @param method DIF detection method; It can be \code{"wald"} for Hou, de la Torre, and Nandakumar's (2014)
#' Wald test method, and \code{"LR"} for likelihood ratio test (Ma, Terzi, Lee,& de la Torre, 2017).
#' @param p.adjust.methods adjusted p-values for multiple hypothesis tests. This is conducted using \code{p.adjust} function in \pkg{stats},
#'  and therefore all adjustment methods supported by \code{p.adjust} can be used, including \code{"holm"},
#'  \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"} and \code{"BY"}. See \code{p.adjust}
#'  for more details. \code{"bonferroni"} is the default.
#' @param anchor.items which items are used as anchors? Default is \code{"none"}, which means none of the items are used as anchors.
#'  For LR method, it can also be an integer vector giving the item numbers for anchors or \code{"all"}, which means all items are treated as anchor items.
#' @param dif.items Items for the DIF detection. By default, all items will be examined.
#' @param approx Whether an approximated LR test is implemented? If TRUE, parameters of items except the studied one will not be re-estimated.
#' @param SE.type Type of standard error estimation methods for the Wald test.
#' @param ... arguments passed to GDINA function for model calibration
#' @return A data frame giving the Wald statistics and associated p-values.
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @seealso \code{\link{GDINA}}
#' @export
#' @examples
#' \dontrun{
#' set.seed(123456)
#' N <- 3000
#' Q <- sim10GDINA$simQ
#' gs <- matrix(c(0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2),ncol = 2, byrow = TRUE)
#' # By default, individuals are simulated from uniform distribution
#' # and deltas are simulated randomly
#' sim1 <- simGDINA(N,Q,gs.parm = gs,model="DINA")
#' sim2 <- simGDINA(N,Q,gs.parm = gs,model=c(rep("DINA",9),"DINO"))
#' dat <- rbind(extract(sim1,"dat"),extract(sim2,"dat"))
#' gr <- c(rep(1,N),rep(2,N))
#' dif.out <- dif(dat,Q,group=gr)
#' dif.out2 <- dif(dat,Q,group=gr,method="LR")
#'}
#' @references
#' Hou, L., de la Torre, J., & Nandakumar, R. (2014). Differential item functioning assessment in cognitive diagnostic modeling: Application of the Wald test to
#' investigate DIF in the DINA model. \emph{Journal of Educational Measurement, 51}, 98-125.
#'
#' Ma, W., Terzi, R., Lee, S., & de la Torre, J. (2017, April). Multiple group cognitive diagnosis models and their applications in detecting differential item functioning. Paper presented at the Annual Meeting ofthe American Educational Research Association, San Antonio, TX.
#'


dif <- function(dat, Q, group, method = "wald", anchor.items = "none", dif.items = "all", p.adjust.methods = "bonferroni", approx = FALSE,
                SE.type = 2, ...){

  if (!is.matrix(dat))
    dat <- as.matrix(dat)
  # default.p.args <- list(initial.anchor.items = "none", alpha.level = .05, maxit = 10, verbose = FALSE)
  # purification.args <- modifyList(default.p.args,purification.args)
  rownames(dat) <- colnames(dat) <- NULL
  est <- NULL
  if (nrow(dat) != length(group))
    stop("The length of group variable must be equal to the number of individuals.",
         call. = FALSE)
  if (length(unique(group)) != 2)
    stop("Only two group DIF can be examined.", call. = FALSE)
  gr <- group
  gr.label <- unique(gr)
  J <- nrow(Q)
  if (length(dif.items) == 1 && dif.items == "all") {
    dif.items <- 1:J
  } else{
    if (min(dif.items) <= 0 |
        max(dif.items) > J)
      stop("dif.items needs to be correctly specified.", call. = FALSE)
    dif.items <- dif.items
  }
  purification.args <- dots("purification.args",list(initial.anchor.items = "none", alpha.level = .05, maxit = 10, verbose = FALSE),...)
  log.purification <- NULL
  # est: unique item parameters for each group on each item
  est <- GDINA::GDINA(dat = bdiagMatrix(list(dat[which(gr==gr.label[1]),],dat[which(gr==gr.label[2]),]),NA),
                      Q = rbind(Q,Q), group = gr,verbose = 0, ...)

  if(tolower(method)=="wald"){

    est1 <- GDINA::GDINA(dat = dat[which(gr==gr.label[1]),], Q = Q,
                         catprob.parm = est$catprob.parm[1:J],
                         att.prior = est$posterior.prob[1,],att.dist = "fixed",verbose = 0,control=list(maxitr = 0))
    est2 <- GDINA::GDINA(dat = dat[which(gr==gr.label[2]),], Q = Q,
                         catprob.parm = est$catprob.parm[(J+1):(2*J)],
                         att.prior = est$posterior.prob[2,],att.dist = "fixed",verbose = 0,control=list(maxitr = 0))

    output <- matrix(0,length(dif.items),3)


      dcov1 <- extract(est1,"delta.cov",SE.type=SE.type)
      dcov2 <- extract(est2,"delta.cov",SE.type=SE.type)
      for (j in 1:length(dif.items)){
        x <- c(extract(est1,"delta.parm")[[dif.items[j]]],
               extract(est2,"delta.parm")[[dif.items[j]]])
        R <- cbind(diag(length(x)/2),-1*diag(length(x)/2))
        vcov <- bdiagMatrix(list(dcov1$cov[dcov1$index$loc[dcov1$index$item==dif.items[j]],
                                    dcov1$index$loc[dcov1$index$item==dif.items[j]]],
                           dcov2$cov[dcov2$index$loc[dcov2$index$item==dif.items[j]],
                                    dcov2$index$loc[dcov2$index$item==dif.items[j]]]))
        output[j,1] <- t(R%*%x)%*%MASS::ginv(R%*%vcov%*%t(R))%*%(R%*%x)
        output[j,2] <- nrow(R)
        output[j,3] <- pchisq(output[j,1],nrow(R),lower.tail = FALSE)
      }


    output <- data.frame(output)
    colnames(output) <- c("Wald stat.","df","p.value")
    rownames(output) <- extract(est,"item.names")[dif.items]
    # output <- Wp
  }else if(method=="LR"){

    if(length(anchor.items)==1 && tolower(anchor.items) == "auto") {
      if(length(purification.args$initial.anchor.items)==1 && tolower(purification.args$initial.anchor.items)=="none"){
        anchoritems.loc <- 0
      }else if(length(purification.args$initial.anchor.items)==1 && tolower(purification.args$initial.anchor.items)=="all"){
        anchoritems.loc <- seq_len(J)
      }else{
        anchoritems.loc <- purification.args$initial.anchor.items
      }
      output <- NULL
      it <- 0

      while(it<purification.args$maxit){
        output <- LRDIF(dat,Q,gr,anchor.items = anchoritems.loc,LR.approx = approx,...)
        it <- it + 1
        log.purification[[it]] <- output
        if(purification.args$verbose){
          cat("Iter = ", it,"Anchor items = Items ",anchoritems.loc,"\n")
          print(output)
        }
        new.anchoritems.loc <- which(output$p.value > purification.args$alpha.level)
        if(length(new.anchoritems.loc)==0)
          new.anchoritems.loc <- 0
        if(identical(anchoritems.loc,new.anchoritems.loc))
          break

        anchoritems.loc <- new.anchoritems.loc
      }
    }else{
      output <- LRDIF(dat,Q,gr,anchor.items = anchor.items,LR.approx = approx,...)
    }

    rownames(output) <- extract(est,"item.names")[dif.items]
    # output <- lr.out
  }
  output$'adj.pvalue' <- stats::p.adjust(output$'p.value', method = p.adjust.methods)
output <- list(test=output,group=gr,p.adjust.methods=p.adjust.methods,TGest = est,log.purification = log.purification)
class(output) <- "dif"
invisible(output)

}



LRDIF <- function(dat, Q, group, anchor.items = "none", dif.items = "all", LR.approx = FALSE,...){

  J <- ncol(dat)
  JJ <- seq_len(J)

  est <- NULL

  gr <- group
  gr.label <- unique(gr)
  J <- nrow(Q)
  if (length(dif.items) == 1){
    if(tolower(dif.items) == "all") {
      dif.items <- JJ
    } else if (min(dif.items) <= 0 |
          max(dif.items) > J){
      stop("difitem needs to be correctly specified.", call. = FALSE)
    }
  }

  if(length(anchor.items)>1 && length(setdiff(JJ, anchor.items))==0) anchor.items <- "all"

  maxit <- 2000
  output <- data.frame(neg2LL=rep(NA,length(dif.items)),
                       LRstat=rep(NA,length(dif.items)),
                       df=rep(NA,length(dif.items)),
                       'p.value'=rep(NA,length(dif.items)))


  if(length(anchor.items)==1 && tolower(anchor.items) == "all") {
    # dif item has dif pars; all other items have same pars for two gr
    if(LR.approx)
      maxit <- c(rep(0,J-1),2000,2000)

    #est: all items have the same pars across groups
    est <- GDINA::GDINA(dat, Q, group = gr,verbose = 0,...)
    item.parm <- extract(est,"catprob.parm")
    for (j in seq_len(length(dif.items))){

      locj <- setdiff(JJ,dif.items[j]) #item no. except the studied one
      estj <- GDINA::GDINA(dat = cbind(dat[,locj],
                                       bdiagMatrix(list(dat[gr==gr.label[1],dif.items[j]],
                                                  dat[gr==gr.label[2],dif.items[j]]),NA)),
                           Q = Q[c(locj,dif.items[j],dif.items[j]),],
                           group = gr,control=list(maxitr = maxit),verbose = 0,
                           catprob.parm = item.parm[c(locj,dif.items[j],dif.items[j])],
                           att.prior =  t(extract(est,"posterior.prob")),...)
      output$LRstat[j] <-  deviance(est) - deviance(estj)
      output$df[j] <- npar(estj)$`No. of total item parameters` - npar(est)$`No. of total item parameters` # distribution parameters identical
      output$neg2LL[j] <- deviance(estj)
      output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
    }
  }else if(length(anchor.items)==1 && (tolower(anchor.items) == "none" || anchor.items == 0)){#dif item has the same par; all other items have dif pars
    # est: unique item parameters for each group on each item
    est <- GDINA::GDINA(dat = bdiagMatrix(list(dat[which(gr==gr.label[1]),],
                                               dat[which(gr==gr.label[2]),]),NA),
                        Q = rbind(Q,Q), group = gr,verbose = 0, ...)
    item.parm <- extract(est,"catprob.parm")
    if(LR.approx)
      maxit <- c(rep(0,2*J-2),2000)


    for (j in seq_len(length(dif.items))){
      locj <- setdiff(JJ,dif.items[j]) #item no. except the studied one
      estj <- GDINA::GDINA(dat = cbind(bdiagMatrix(list(dat[gr==gr.label[1],locj],
                                                  dat[gr==gr.label[2],locj]),NA),
                                       dat[,dif.items[j]]),
                           Q = Q[c(locj,locj,dif.items[j]),],
                           group = gr,control=list(maxitr = maxit),verbose = 0,
                           catprob.parm = item.parm[c(locj,(locj+J),dif.items[j])],
                           att.prior = t(extract(est,"posterior.prob")),...)
      output$LRstat[j] <- deviance(estj) - deviance(est)
      output$df[j] <- npar(est)$`No. of total item parameters` - npar(estj)$`No. of total item parameters`
      output$neg2LL[j] <- deviance(estj)
      output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
    }


  }else{
    variant.items <- setdiff(JJ, anchor.items)
    # est: anchor items + non-anchor for g1 + non-anchor for g2
    est <- GDINA::GDINA(dat = cbind(dat[,anchor.items],
                                    bdiagMatrix(list(dat[which(gr==gr.label[1]),variant.items],
                                                     dat[which(gr==gr.label[2]),variant.items]),
                                                NA)),
                        Q = Q[c(anchor.items,variant.items,variant.items),],
                        group = gr,verbose = 0,...)
    item.parm <- extract(est,"catprob.parm")
    item.set <- c(anchor.items,variant.items,variant.items)


    for (j in seq_len(length(dif.items))){
      # if the studied item is an anchor item, it will be treated as a non-anchor item
      if(dif.items[j]%in%anchor.items){
        if(LR.approx)
          maxit <- c(rep(0,length(item.set)-1),rep(2000,2))

        locj <- which(dif.items[j]==item.set)
        updated.anchor.items <- anchor.items[-locj]
        estj <- GDINA::GDINA(dat = cbind(dat[,updated.anchor.items],
                                        bdiagMatrix(list(dat[which(gr==gr.label[1]),variant.items],
                                                         dat[which(gr==gr.label[2]),variant.items]),
                                                    NA),
                                        bdiagMatrix(list(dat[which(gr==gr.label[1]),dif.items[j]],
                                                         dat[which(gr==gr.label[2]),dif.items[j]]),
                                                    NA)),
                            Q = Q[c(updated.anchor.items,variant.items,variant.items,
                                    dif.items[j],dif.items[j]),],
                            group = gr,verbose = 0, control=list(maxitr = maxit),
                            catprob.parm = item.parm[c(seq_len(length(item.parm))[-locj],locj,locj)],
                            att.prior = t(extract(est,"posterior.prob")),...)
        output$LRstat[j] <- deviance(est) - deviance(estj)
        output$df[j] <- npar(estj)$`No. of total item parameters` - npar(est)$`No. of total item parameters`
        output$neg2LL[j] <- deviance(estj)
        output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
      }else{

        # if the studied item is not an anchor item, it will be treated as an anchor item
        if(LR.approx)
          maxit <- c(rep(0,length(item.set)-2),2000)

        if(length(variant.items)>1){

          locj <- which(variant.items!=dif.items[j])

          estj <- GDINA::GDINA(dat = cbind(dat[,anchor.items],
                                           bdiagMatrix(list(dat[which(gr==gr.label[1]),variant.items[locj]],
                                                            dat[which(gr==gr.label[2]),variant.items[locj]]),
                                                       NA),
                                           dat[,dif.items[j]]),
                               Q = Q[c(anchor.items,variant.items[locj],variant.items[locj],dif.items[j]),],
                               group = gr,verbose = 0, control=list(maxitr = maxit),
                               catprob.parm = item.parm[c(which(item.set!=dif.items[j]),which(item.set==dif.items[j])[1])],
                               att.prior = t(extract(est,"posterior.prob")),...)
          output$LRstat[j] <- deviance(estj) - deviance(est)
          output$df[j] <- npar(est)$`No. of total item parameters` - npar(estj)$`No. of total item parameters`
          output$neg2LL[j] <- deviance(estj)
          output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
        }else{ # only one variant item

          estj <- GDINA::GDINA(dat = dat[,c(anchor.items,variant.items)],
                               Q = Q[c(anchor.items,variant.items),],
                               group = gr,verbose = 0, control=list(maxitr = maxit),
                               catprob.parm = item.parm[c(which(item.set!=dif.items[j]),which(item.set==dif.items[j])[1])],
                               att.prior = t(extract(est,"posterior.prob")),...)
          output$LRstat[j] <- deviance(estj) - deviance(est)
          output$df[j] <- npar(est)$`No. of total item parameters` - npar(estj)$`No. of total item parameters`
          output$neg2LL[j] <- deviance(estj)
          output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
        }

      }

  }


  }
  output
}

