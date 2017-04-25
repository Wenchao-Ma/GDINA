#' @include GDINA.R
#' @title  Differential item functioning for cognitive diagnosis models
#'
#' @description   This function is used to detect differential item functioning based on the models estimated
#' in the \code{\link{GDINA}} function using the Wald test (Hou, de la Torre, & Nandakumar, 2014)
#' @param method DIF detection method; It can be \code{"wald"} for Hou, de la Torre, and Nandakumar's (2014)
#' Wald test method, and \code{"LR"} for likelihood ratio test.
#' @param parm The type of parameters associated with the Wald test for the DIF detection. It can be either \code{"itemprob"}
#'  or \code{"delta"} for item probabilities and delta parameters, respectively.
#' @param p.adjust.methods adjusted p-values for multiple hypothesis tests. This is conducted using \code{p.adjust} function in \pkg{stats},
#'  and therefore all adjustment methods supported by \code{p.adjust} can be used, including \code{"holm"},
#'  \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"} and \code{"BY"}. See \code{p.adjust}
#'  for more details. \code{"bonferroni"} is the default.
#' @param difitem Items for the DIF detection. By default, all items will be examined.
#' @param LR.type Type of likelihood ratio test for DIF detection. It can be \code{'free.all'} or
#' \code{'free.one'}.
#' @param SE.type Type of standard error estimation methods for the Wald test.
#' @inheritParams GDINA
#' @param ... Other arguments passed to GDINA function for model calibration
#' @return A data frame giving the Wald statistics and associated p-values.
#'
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @seealso \code{\link{GDINA}}
#' @export
#' @examples
#' \dontrun{
#' set.seed(123456)
#' N <- 1000
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
#' gr <- c(rep("male",N),rep("female",N))
#' dif.out <- dif(dat,Q,group=gr)
#' plotIRF(dif.out,4)
#' dif.out2 <- dif(dat,Q,group=gr,method="LR")
#'}
#' @references
#' Hou, L., de la Torre, J., & Nandakumar, R. (2014). Differential item functioning assessment in cognitive diagnostic modeling: Application of the Wald test to
#' investigate DIF in the DINA model. \emph{Journal of Educational Measurement, 51}, 98-125.
#'


dif <- function(dat, Q, group, method = "wald", p.adjust.methods = "bonferroni",LR.type="free.all",
                difitem = "all", parm = "delta", digits = 4, SE.type = 2,...){
  if (length(group)==1){
    gr <- dat[,group]
    dat <- dat[,-group]
  }else{
    if (nrow(dat)!=length(group))stop("The length of group variable must be equal to the number of individuals.",call. = FALSE)
  if (length(unique(group))!=2)stop("Only two group DIF can be examined.",call. = FALSE)
  gr <- group
  }
  gr.label <- unique(gr)
  J <- nrow(Q)
  if(length(difitem)==1&&difitem == "all") difitem <- 1:J
  est <- GDINA::GDINA(dat = bdiagMatrix(list(dat[which(gr==gr.label[1]),],dat[which(gr==gr.label[2]),]),NA),
                      Q = rbind(Q,Q),
                      group = gr,...)

  if(method=="wald"){

    output <- matrix(0,length(difitem),3)

    if(parm=="itemprob"){
      pcov <- extract(est,"catprob.cov",SE.type=SE.type)
      for (j in 1:length(difitem)){
        # category probabilities
        if(extract(est,"models_numeric")[[difitem[j]]]%in%c(1,2)){
          x <- c(itemparm(est,"gs",digits = 8)[difitem[j],],
                 itemparm(est,"gs",digits = 8)[difitem[j]+J,])
        }else{
          x <- c(extract(est,"catprob.parm")[[difitem[j]]],
                 extract(est,"catprob.parm")[[difitem[j]+J]])
        }

        R <- cbind(diag(length(x)/2),-1*diag(length(x)/2))
        vcov <- bdiagMatrix(list(pcov$cov[pcov$index$loc[pcov$index$item==difitem[j]],
                                    pcov$index$loc[pcov$index$item==difitem[j]]],
                           pcov$cov[pcov$index$loc[pcov$index$item==difitem[j]+J],
                                    pcov$index$loc[pcov$index$item==difitem[j]+J]]))
        output[j,1] <- t(R%*%x)%*%MASS::ginv(R%*%vcov%*%t(R))%*%(R%*%x)
        output[j,2] <- nrow(R)
        output[j,3] <- pchisq(output[j,1],nrow(R),lower.tail = FALSE)
      }
    }else if(parm=="delta"){
      dcov <- extract(est,"delta.cov",SE.type=SE.type)
      for (j in 1:length(difitem)){
        x <- c(extract(est,"delta.parm")[[difitem[j]]],
               extract(est,"delta.parm")[[difitem[j]+J]])
        R <- cbind(diag(length(x)/2),-1*diag(length(x)/2))
        vcov <- bdiagMatrix(list(dcov$cov[dcov$index$loc[dcov$index$item==difitem[j]],
                                    dcov$index$loc[dcov$index$item==difitem[j]]],
                           dcov$cov[dcov$index$loc[dcov$index$item==difitem[j]+J],
                                    dcov$index$loc[dcov$index$item==difitem[j]+J]]))
        output[j,1] <- t(R%*%x)%*%MASS::ginv(R%*%vcov%*%t(R))%*%(R%*%x)
        output[j,2] <- nrow(R)
        output[j,3] <- pchisq(output[j,1],nrow(R),lower.tail = FALSE)
      }

    }
    output <- round(data.frame(output),digits)
    colnames(output) <- c("Wald stat.","df","p.value")
    rownames(output) <- extract(est,"item.names")[difitem]
    # output <- Wp
  }else if(method=="LR"){
    output <- data.frame(neg2LL=rep(NA,length(difitem)),
                         LRstat=rep(NA,length(difitem)),
                         df=rep(NA,length(difitem)),
                         'p.value'=rep(NA,length(difitem)))

     if(LR.type=="free.one") {
       # dif item has dif pars; all other items have same pars for two gr
       est0 <- GDINA::GDINA(dat, Q, group = gr,...)

        for (j in 1:length(difitem)){

            gr1dat <- gr2dat <- dat[,difitem[j]]
            gr1dat[gr==gr.label[2]] <- NA
            gr2dat[gr==gr.label[1]] <- NA
            estj <- GDINA::GDINA(dat = cbind(gr1dat,gr2dat,dat[,-difitem[j]]),
                                 Q = rbind(Q[difitem[j],],Q[difitem[j],],Q[-difitem[j],]),
                                 group = gr,
                                 catprob.parm = c(list(extract(est0,"catprob.parm")[[difitem[j]]],
                                                       extract(est0,"catprob.parm")[[difitem[j]]]),
                                                  extract(est0,"catprob.parm")[-difitem[j]]),
                                 att.prior = extract(est0,"att.prior"),...)
            output$LRstat[j] <-  deviance(est0) - deviance(estj)
            output$df[j] <- npar(estj)$`No. of parameters` - npar(est0)$`No. of parameters`
            output$neg2LL[j] <- deviance(estj)
            output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
        }
      }else if(LR.type=="free.all"){#dif item has the same par; all other items have dif pars
        # gr1dat <- gr2dat <- dat[,-difitem[j]]
        # gr1dat[gr==gr.label[2],] <- NA
        # gr2dat[gr==gr.label[1],] <- NA
        # newQ <- rbind(Q[difitem[j],],Q[-difitem[j],],Q[-difitem[j],])
        # newdat <- cbind(dat[,difitem[j]],gr1dat,gr2dat)
        # tmp.par <- c(list(extract(est,"catprob.parm")[[difitem[j]]]),
        #              extract(est,"catprob.parm")[-c(difitem[j],difitem[j]+J)])
        for (j in 1:length(difitem)){
          estj <- GDINA::GDINA(dat = cbind(dat[,difitem[j]],bdiag(list(dat[gr==gr.label[1],-difitem[j]],dat[gr==gr.label[2],-difitem[j]]),NA)),
                             Q = rbind(Q[difitem[j],],Q[-difitem[j],],Q[-difitem[j],]),
                             group = gr,
                             catprob.parm = c(list(extract(est,"catprob.parm")[[difitem[j]]]),
                                              extract(est,"catprob.parm")[-c(difitem[j],difitem[j]+J)]),
                             att.prior = extract(est,"att.prior"),...)
        output$LRstat[j] <- deviance(estj) - deviance(est)
        output$df[j] <- npar(est)$`No. of parameters` - npar(estj)$`No. of parameters`
        output$neg2LL[j] <- deviance(estj)
        output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
      }


    }
    output <- round(output,digits)
    rownames(output) <- extract(est,"item.names")[difitem]
    # output <- lr.out
  }
  output$'adj.pvalue' <- stats::p.adjust(output$'p.value', method = p.adjust.methods)
output <- list(test=output,group=gr,mg.est=est,p.adjust.methods=p.adjust.methods)
class(output) <- "dif"
return(output)

}
