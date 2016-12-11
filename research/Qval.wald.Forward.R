#' Qvalidation for sequential G-DINA model using forward search (under development)
#'
#' @param GDINA.obj GDINA object
#' @param item items for Q-matrix validation
#' @param SE.type Typeof standard errors
#' @param alpha.level alpha level
#' @param GDI GDI or not
#' @param ... other arguments
#'
#' @return a list
# #' @export
#'
#' @examples
#'
#' est <- GDINA(sim21seqDINA$simdat,sim21seqDINA$simQc,sequential = TRUE)
#' x <- Qval.wald.forward(est)
#'
Qval.wald.forward <- function(GDINA.obj, item = NULL, SE.type = 2, alpha.level = 0.05, GDI = TRUE, ...){

  dat <- internalextract(GDINA.obj,"dat")

  Qr <- Qc <- internalextract(GDINA.obj,"Qc")

  Q <- internalextract(GDINA.obj,"Q")

  seqdat <- seq_coding(dat, Qc)

  N <- internalextract(GDINA.obj,"nobs")

  # number of categories
  J <- internalextract(GDINA.obj,"ncat")

  K <- internalextract(GDINA.obj,"natt")

  L <- no_LC(Q)

  Kj <- rowSums(alpha(K)[-1,])
  w <- internalextract(GDINA.obj,"posterior.prob") #1 x L

  RN <- NgRg(GDINA.obj$technicals$logposterior.i,seqdat,eta.loc(matrix(1,J,K)),1-is.na(seqdat))
  expectedR <- RN$Rg
  expectedN <- RN$Ng
  # rc <- apply(seqdat,2,function(x){
  #   colSums(x*exp(internalextract(GDINA.obj,"logposterior.i")),na.rm = TRUE)
  # })

  # est.p <- rc/(w*N)
  est.p <- t((expectedR)/(expectedN))
  loc <- eta.loc(diag(K)) #2^K-1 x 2^K
   vsg <- varsigma(as.matrix(t(loc)),as.matrix(est.p),c(w))

  #first attribute is the one with the largest GDI
  first.att <- apply(vsg,1,which.max.randomtie)
  att.monitor <- Qrr <- vector("list",J)
  inichoose <- numeric(J)
  fullset <- c(1:K)
  # iteras <- NULL
  if (is.null(item)) item <- c(1:nrow(Q))
  for (j in item) {

    inichoose[j] <- currentset <- first.att[j] # initial att. ---largest GDI
    att.monitor[[j]] <- c(att.monitor[[j]],currentset)
      loop <- T
      it <- 1
      while(loop&&it<K){
        difset <- setdiff(fullset,currentset)
        add.a <- NULL
        #********************************************************************Second round eval. Wald forward
        Rm <- Rmatrix.att(length(currentset)+1)
        Wp.a <- NULL
        for (k in difset){ # 2nd att.
          cat("Item",j,"Att",k,"\n")
          Qcr <- Qc
          Qcr[j,3:ncol(Qcr)] <- 0
          Qcr[j,2+c(currentset,k)] <- 1
          etaj <- eta.loc(Qcr[,-c(1:2)])[j,]
          itemparj <- internalextract(GDINA.obj,"itemprob.parm")
          itemparj[[j]] <- aggregate(expectedR[j,],by=list(etaj),sum)$x/
            aggregate(expectedN[j,],by=list(etaj),sum)$x
          out <- GDINA(dat,Qcr,maxitr = 0, sequential = TRUE,
                       itemprob.parm = itemparj,
                       att.prior = exp(internalextract(GDINA.obj,"posterior.prob")))
          cov <- internalextract(out,"itemprob.cov",type=SE.type)
          covind <- cov$index
          Varj <- cov$cov[covind[which(covind$item==j),"loc"],covind[which(covind$item==j),"loc"]]
          #print(V$SE[[j]])
          # item pars
          Cparj <- itemparm(out)[[j]]

         Wp.a <- rbind(Wp.a,sapply(Rm,function(x) {
            # wald statistic & p values
            W <-t(x %*% Cparj) %*% MASS::ginv(x %*% Varj %*% t(x)) %*% (x %*% Cparj)
            p <- pchisq(W,nrow(x),lower.tail = F)
            return (p)
          }))

        }
        #********************************************************************end 2nd round eval. att. Forward
        add.new <- NULL
        for (d in 1:length(difset)){
          add.new <- rbind(add.new,c(difset[d],Wp.a[d,1+sum(difset[d]>currentset)],vsg[j,difset[d]]))
        }
        # iteras <- rbind(iteras,cbind(j,it,add.new))
        #*********************************************************If additional elements should be added
        if (any(add.new[,2]<alpha.level)){ # Yes
          add.new <- add.new[which(add.new[,2]<alpha.level),,drop=FALSE]

          if (GDI) add.a <- add.new[which.max(add.new[,3]),1] # the one with largest GDI
          if (!GDI) add.a <- add.new[which.max(add.new[,2]),1] # the one with largest p-value


          currentset <- sort(c(currentset,add.a)) # all att. chosen
          att.monitor[[j]] <- c(att.monitor[[j]],add.a)
          loop <- T
          it <- it + 1
        }else{   # No
          loop <- F
        }
      }

      Qr[j,3:ncol(Qr)] <- 0
      Qr[j,2+currentset] <- 1


    }
    #print(do.call(rbind,Qrr))

  return(list(sugQ=Qr,initialAtt=first.att,Qc=Qc))
  }





