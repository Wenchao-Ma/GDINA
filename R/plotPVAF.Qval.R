#' Mesa plot for Q-matrix validation
#'
#' @param Qval.obj model object of class \code{Qvalidation}
#' @param item a vector specifying which item(s) the plots are drawn for
#' @param type types of the plot. It can be \code{"best"} or \code{"all"}. If \code{"best"},
#'     for all q-vectors requiring the same number of attributes, only the one with the largest PVAF
#'     is plotted, which means \eqn{K_j} q-vectors are plotted; If \code{"all"}, all q-vectors
#'     will be plotted.
#' @param no.qvector the number of q vectors that need to be plotted when \code{type="all"}. The default is 10,
#'        which means the 10 q vectors with the largest PVAFs are plotted.
#' @param data.label logical; To show data label or not?
#' @param original.q.label logical; print the label showing the original q-vector or not?
#' @param auto.ylim logical; create y range automatically or not?
#' @param ... additional arguments
#' @seealso \code{\link{Qval}}, \code{\link{autoGDINA}}
#' @examples
#'
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' Q[1,] <- c(0,1,0)
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' out <- Qval(mod1,eps = 0.9)
#' item <- c(1,2,10)
#' mesaplot(out,item=item,data.label=TRUE,type="all")
#' mesaplot(out,item=10,type="best")
#' mesaplot(out,item=10,type="all",no.qvector=5)
#'
#'
#'
#' @export
mesaplot <- function(Qval.obj, item, type = "best", no.qvector = 10,
                     data.label = FALSE,
                     original.q.label = TRUE,auto.ylim = TRUE,...){
  UseMethod("mesaplot")
}
#' @export
mesaplot.Qval <-
  function(Qval.obj, item, type = "best", no.qvector = 10,
           data.label = FALSE,
           original.q.label = TRUE,auto.ylim = TRUE,...)
  {
    Q <- extract.Qval(Qval.obj,"Q")
    # Q <- Qval.obj$sugg.Q[,grep("orig.q",colnames(Qvalidation.obj$sugg.Q))]
    K <- ncol(Q)
    L <- (2^K-1) # L-1
    fullPVAF <- extract.Qval(Qval.obj,"PVAF")
    if(tolower(type)=="all"){

      if (L<no.qvector) no.qvector <- L
      for (y in item){
        #which one is the true
        locy0 <- which(apply(alpha(K)[-1,],1,function(x){
          all(x==Q[y,])}))

        locy <- no.qvector-(L-which(order(fullPVAF[,y],decreasing = F)==locy0))
        ordered.PVAF.j <- sort(fullPVAF[,y],decreasing = FALSE)
        plot(ordered.PVAF.j[(L-no.qvector+1):L],xaxt="n",type="o",ylab = "PVAF",xlab="q-vectors",
             main = paste("Mesa Plot for Item",y),ylim = c(0,1))
        axis(1,at=c(1:no.qvector),labels = names(ordered.PVAF.j[(L-no.qvector+1):L]))
        if (locy>0){
          points(locy,fullPVAF[locy0,y],col="red",pch=19)
        }
        if (original.q.label) text(no.qvector-1,0.15,paste("original q-vector:\n",names(fullPVAF[,y])[locy0]))
        if (auto.ylim) ylim = c(max(0,round(min(ordered.PVAF.j)-0.1,1)),1) else ylim=c(0,1)
        yloc <- ordered.PVAF.j[(L-no.qvector+1):L]-diff(ylim)/15
        yloc[yloc<=ylim[1]] <- yloc[yloc<=ylim[1]] + 2 * diff(ylim)/15
        if (data.label) text(c(1:no.qvector),
                             yloc,
                             ordered.PVAF.j[(L-no.qvector+1):L])

      }
    }else if(tolower(type)=="best"){
      Kj <- rowSums(alpha(K)[-1,])
      bestPVAF <- aggregate(fullPVAF,by=list(Kj),max)[,-1]
      label.bestPVAF <- rownames(fullPVAF)
      for(j in item){
        bestloc <- match(bestPVAF[,j],fullPVAF[,j])
        if (auto.ylim) ylim = c(max(0,round(min(bestPVAF[,j])-0.1,1)),1) else ylim=c(0,1)
        plot(bestPVAF[,j],xaxt="n",type="o",ylab = "PVAF",xlab="q-vectors",
             main = paste("Mesa Plot for Item",j),ylim = ylim)
        axis(1,at=c(1:nrow(bestPVAF)),labels = label.bestPVAF[bestloc])
        yloc <- bestPVAF[,j]-diff(ylim)/15
        yloc[yloc<=ylim[1]] <- yloc[yloc<=ylim[1]] + 2 * diff(ylim)/15
        if (data.label) text(c(1:nrow(bestPVAF)),yloc,bestPVAF[,j])
        locy0 <- which(apply(alpha(K)[-1,],1,function(x){
          all(x==Q[j,])}))
        if(locy0%in%bestloc) points(which(bestloc==locy0),fullPVAF[locy0,j],col="red",pch=19)
        if (original.q.label) text(K-1,ylim[1]+diff(ylim)/6,paste("original q-vector:\n",names(fullPVAF[,j])[locy0]))
      }
    }

  }

