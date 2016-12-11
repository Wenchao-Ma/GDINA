#' Differential item functioning for cognitive diagnosis models
#'
#' This function is used to detect differential item functioning based on the models estimated
#' in the \code{\link{GDINA}} function using the Wald test (Hou, de la Torre, & Nandakumar, 2014)
#' @param method DIF detection method; It can be \code{"wald"} for Hou, de la Torre, and Nandakumar's (2014)
#' Wald test method, and \code{"LR"} for likelihood ratio test.
#' @param parm The type of parameters associated with the Wald test for the DIF detection. It can be either \code{"itemprob"}
#'  or \code{"delta"} for item probabilities and delta parameters, respectively.
#' @param difitem Items for the DIF detection. By default, all items will be examined.
#' @param SE.type Type of standard error estimation methods for Wald test.
#' @inheritParams GDINA
#' @param ... arguments passed to GDINA function for model calibration
#' @return a data frame giving the Wald statistics and associated p-values.
#'
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @seealso \code{\link{GDINA}}
#' @export
#' @examples
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
#'
#' @references
#' Hou, L., de la Torre, J., & Nandakumar, R. (2014). Differential item functioning assessment in cognitive diagnostic modeling: Application of the Wald test to
#' investigate DIF in the DINA model. \emph{Journal of Educational Measurement, 51}, 98-125.
#'


dif <- function(dat, Q, group, method = "wald", difitem = "all", parm = "delta", digits = 4, SE.type = 2,...){
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
  if(difitem == "all") difitem <- 1:J
gr1dat <- gr2dat <- dat
gr1dat[gr==gr.label[2],] <- NA
gr2dat[gr==gr.label[1],] <- NA
est <- GDINA::GDINA(cbind(gr1dat,gr2dat), rbind(Q,Q), group = gr,...)

  if(method=="wald"){

    output <- matrix(0,length(difitem),3)

    if(parm=="itemprob"){
      pcov <- internalextract(est,"catprob.cov",type=SE.type)
      for (j in difitem){
        # category probabilities
        if(internalextract(est,"models_numeric")[[j]]%in%c(1,2)){
          x <- c(itemparm(est,"gs",digits = 8)[j,],itemparm(est,"gs",digits = 8)[j+J,])
        }else{
          x <- c(internalextract(est,"catprob.parm")[[j]],internalextract(est,"catprob.parm")[[j+J]])
        }

        R <- cbind(diag(length(x)/2),-1*diag(length(x)/2))
        vcov <- bdiag(list(pcov$cov[pcov$index$loc[pcov$index$item==j],pcov$index$loc[pcov$index$item==j]],
                           pcov$cov[pcov$index$loc[pcov$index$item==j+J],pcov$index$loc[pcov$index$item==j+J]]))
        output[j,1] <- t(R%*%x)%*%MASS::ginv(R%*%vcov%*%t(R))%*%(R%*%x)
        output[j,2] <- nrow(R)
        output[j,3] <- pchisq(output[j,1],nrow(R),lower.tail = FALSE)
      }
    }else if(parm=="delta"){
      dcov <- internalextract(est,"delta.cov",type=SE.type)
      for (j in difitem){
        x <- c(internalextract(est,"delta.parm")[[j]],internalextract(est,"delta.parm")[[j+J]])
        R <- cbind(diag(length(x)/2),-1*diag(length(x)/2))
        vcov <- bdiag(list(dcov$cov[dcov$index$loc[dcov$index$item==j],dcov$index$loc[dcov$index$item==j]],
                           dcov$cov[dcov$index$loc[dcov$index$item==j+J],dcov$index$loc[dcov$index$item==j+J]]))
        output[j,1] <- t(R%*%x)%*%MASS::ginv(R%*%vcov%*%t(R))%*%(R%*%x)
        output[j,2] <- nrow(R)
        output[j,3] <- pchisq(output[j,1],nrow(R),lower.tail = FALSE)
      }

    }
    output <- round(data.frame(output),digits)
    colnames(output) <- c("Wald stat.","df","p.value")
    rownames(output) <- internalextract(est,"item.names")[difitem]
    # output <- Wp
  }else if(method=="LR"){
    output <- data.frame(neg2LL=rep(NA,nrow(Q)),
                         LRstat=rep(NA,nrow(Q)),
                         df=rep(NA,nrow(Q)),
                         'p.value'=rep(NA,nrow(Q)))

    for (j in difitem){
      gr1dat <- gr2dat <- dat[,-j]
      gr1dat[gr==gr.label[2],] <- NA
      gr2dat[gr==gr.label[1],] <- NA
      newQ <- rbind(Q[j,],Q[-j,],Q[-j,])
      newdat <- cbind(dat[,j],gr1dat,gr2dat)
      estj <- GDINA::GDINA(newdat, newQ, group = gr,...)
      output$neg2LL[j] <- deviance(estj)
      output$LRstat[j] <- deviance(estj) - deviance(est)
      output$df[j] <- npar(est)[1] - npar(estj)[1]
      output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
    }
    output <- round(output,digits)
    rownames(output) <- internalextract(est,"item.names")[1:nrow(Q)]
    # output <- lr.out
  }
output <- list(test=output,group=gr,mg.est=est)
class(output) <- "dif"
return(output)

}
