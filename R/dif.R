#' @include GDINA.R
#' @title  Differential item functioning for cognitive diagnosis models
#'
#' @description   This function is used to detect differential item functioning based on the models estimated
#' in the \code{\link{GDINA}} function using the Wald test (Hou, de la Torre, & Nandakumar, 2014) and the likelihood ratio
#' test (Ma, Terzi, Lee, & de la Torre, 2017). It can only detect DIF for two groups currently.
#'
#' @param method DIF detection method; It can be \code{"wald"} for Hou, de la Torre, and Nandakumar's (2014)
#' Wald test method, and \code{"LR"} for likelihood ratio test (Ma, Terzi, Lee,& de la Torre, 2017).
#' @param p.adjust.methods adjusted p-values for multiple hypothesis tests. This is conducted using \code{p.adjust} function in \pkg{stats},
#'  and therefore all adjustment methods supported by \code{p.adjust} can be used, including \code{"holm"},
#'  \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"} and \code{"BY"}. See \code{p.adjust}
#'  for more details. \code{"bonferroni"} is the default.
#' @param difitem Items for the DIF detection. By default, all items will be examined.
#' @param LR.type Type of likelihood ratio test for DIF detection. It can be \code{'free.all'} or
#' \code{'free.one'}.
#' @param LR.approx Whether an approximated LR test is implemented? If TRUE, anchor item parameters will not be re-estimated but fixed.
#' @param SE.type Type of standard error estimation methods for the Wald test.
#' @inheritParams GDINA
#' @param ... Other arguments passed to GDINA function for model calibration
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


dif <- function(dat, Q, group, method = "wald", p.adjust.methods = "bonferroni",LR.type="free.all", LR.approx = FALSE,
                difitem = "all", digits = 4, SE.type = 2, ...){
  if(!is.matrix(dat)){dat <- as.matrix(dat)}
  rownames(dat) <- colnames(dat) <- NULL

  if (nrow(dat)!=length(group))stop("The length of group variable must be equal to the number of individuals.",call. = FALSE)
  if (length(unique(group))!=2)stop("Only two group DIF can be examined.",call. = FALSE)
  gr <- group
  gr.label <- unique(gr)
  J <- nrow(Q)
  if(length(difitem)==1&&difitem == "all") {
    difitems <- 1:J
  }else{
    if(min(difitem)<=0|max(difitem)>J) stop("difitem needs to be correctly specified.",call. = FALSE)
    difitems <- difitem
  }


  if(tolower(method)=="wald"){
    # est: unique item parameters for each group on each item
    est <- GDINA::GDINA(dat = bdiagMatrix(list(dat[which(gr==gr.label[1]),],dat[which(gr==gr.label[2]),]),NA),
                        Q = rbind(Q,Q), group = gr,verbose = 0, ...)

    est1 <- GDINA::GDINA(dat = dat[which(gr==gr.label[1]),], Q = Q,
                         catprob.parm = est$catprob.parm[1:J],
                         att.prior = est$posterior.prob[1,],att.dist = "fixed",verbose = 0,control=list(maxitr = 0),...)
    est2 <- GDINA::GDINA(dat = dat[which(gr==gr.label[2]),], Q = Q,
                         catprob.parm = est$catprob.parm[(J+1):(2*J)],
                         att.prior = est$posterior.prob[2,],att.dist = "fixed",verbose = 0,control=list(maxitr = 0),...)

    output <- matrix(0,length(difitems),3)


      dcov1 <- extract(est1,"delta.cov",SE.type=SE.type)
      dcov2 <- extract(est2,"delta.cov",SE.type=SE.type)
      for (j in 1:length(difitems)){
        x <- c(extract(est1,"delta.parm")[[difitems[j]]],
               extract(est2,"delta.parm")[[difitems[j]]])
        R <- cbind(diag(length(x)/2),-1*diag(length(x)/2))
        vcov <- bdiagMatrix(list(dcov1$cov[dcov1$index$loc[dcov1$index$item==difitems[j]],
                                    dcov1$index$loc[dcov1$index$item==difitems[j]]],
                           dcov2$cov[dcov2$index$loc[dcov2$index$item==difitems[j]],
                                    dcov2$index$loc[dcov2$index$item==difitems[j]]]))
        output[j,1] <- t(R%*%x)%*%MASS::ginv(R%*%vcov%*%t(R))%*%(R%*%x)
        output[j,2] <- nrow(R)
        output[j,3] <- pchisq(output[j,1],nrow(R),lower.tail = FALSE)
      }


    output <- round(data.frame(output),digits)
    colnames(output) <- c("Wald stat.","df","p.value")
    rownames(output) <- extract(est,"item.names")[difitems]
    # output <- Wp
  }else if(method=="LR"){
    output <- data.frame(neg2LL=rep(NA,length(difitems)),
                         LRstat=rep(NA,length(difitems)),
                         df=rep(NA,length(difitems)),
                         'p.value'=rep(NA,length(difitems)))

     if(LR.type=="free.one") {
       # dif item has dif pars; all other items have same pars for two gr
       if(LR.approx){
         maxit <- c(1000,1000,rep(0,J-1))
       }else{
         maxit <- 1000
       }
       #est0: all items have the same pars across groups
       est <- est0 <- GDINA::GDINA(dat, Q, group = gr,verbose = 0,...)

        for (j in 1:length(difitems)){

            gr1dat <- gr2dat <- dat[,difitems[j]]
            gr1dat[gr==gr.label[2]] <- NA
            gr2dat[gr==gr.label[1]] <- NA
            estj <- GDINA::GDINA(dat = cbind(gr1dat,gr2dat,dat[,-difitems[j]]),
                                 Q = rbind(Q[difitems[j],],Q[difitems[j],],Q[-difitems[j],]),
                                 group = gr,control=list(maxitr = maxit),verbose = 0,
                                 catprob.parm = c(list(extract(est0,"catprob.parm")[[difitems[j]]],
                                                       extract(est0,"catprob.parm")[[difitems[j]]]),
                                                  extract(est0,"catprob.parm")[-difitems[j]]),
                                 att.prior =  t(extract(est0,"posterior.prob")),...)
            output$LRstat[j] <-  deviance(est0) - deviance(estj)
            output$df[j] <- npar(estj)$`No. of parameters` - npar(est0)$`No. of parameters`
            output$neg2LL[j] <- deviance(estj)
            output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
        }
      }else if(LR.type=="free.all"){#dif item has the same par; all other items have dif pars
        # est: unique item parameters for each group on each item
        est <- GDINA::GDINA(dat = bdiagMatrix(list(dat[which(gr==gr.label[1]),],dat[which(gr==gr.label[2]),]),NA),
                            Q = rbind(Q,Q), group = gr,verbose = 0, ...)
        if(LR.approx){
          maxit <- c(1000,rep(0,2*J-2))
        }else{
          maxit <- 1000
        }

        for (j in 1:length(difitems)){
          estj <- GDINA::GDINA(dat = cbind(dat[,difitems[j]],
                                           bdiag(list(dat[gr==gr.label[1],-difitems[j]],
                                                      dat[gr==gr.label[2],-difitems[j]]),NA)),
                             Q = rbind(Q[difitems[j],],Q[-difitems[j],],Q[-difitems[j],]),
                             group = gr,control=list(maxitr = maxit),verbose = 0,
                             catprob.parm = c(list(extract(est,"catprob.parm")[[difitems[j]]]),
                                              extract(est,"catprob.parm")[-c(difitems[j],difitems[j]+J)]),
                             att.prior = t(extract(est,"posterior.prob")),...)
        output$LRstat[j] <- deviance(estj) - deviance(est)
        output$df[j] <- npar(est)$`No. of parameters` - npar(estj)$`No. of parameters`
        output$neg2LL[j] <- deviance(estj)
        output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
      }


    }
    output <- round(output,digits)
    rownames(output) <- extract(est,"item.names")[difitems]
    # output <- lr.out
  }
  output$'adj.pvalue' <- stats::p.adjust(output$'p.value', method = p.adjust.methods)
output <- list(test=output,group=gr,p.adjust.methods=p.adjust.methods)
class(output) <- "dif"
return(output)

}
