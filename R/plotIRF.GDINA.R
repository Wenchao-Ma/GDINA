#' @include GDINA.R
#' @title  Create plots for GDINA estimates
#'
#' @description   Create various plots for GDINA estimates
#'
#' @param x model object of class \code{\link{GDINA}}
#' @param what type of plot. Can be \code{"IRF"} for item/category response function plot,
#'      or \code{"mp"} for mastery probabilities for individuals.
#' @param IRF.args arguments for IRF plot - a list with the following possible options:
#'  \itemize{
#'    \item \code{item} - A scalar or vector specifying the item(s) for IRF plots.
#'    \item \code{erorbar} - Add error bar (estimate - SE, estimate + SE) to the IRF plots?
#'    \item \code{SE.type} - How is SE estimated. By default, it's based on OPG using incomplete information.
#'    }
#' @param mp.args arguments for mastery probability plot - a list with the following possible options:
#'  \itemize{
#'    \item \code{id} - A scalar or vector specifying the number of individuals for mastery plots.
#'    \item \code{att.names} - Optional; a vector for attribute names.#'    }
#' @param ... additional arguments
#' @seealso \code{\link{GDINA}}, \code{\link{autoGDINA}}
#' @export
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' #plot item response functions for item 10
#' plot(mod1, what = "IRF", IRF.args=list(item = 10))
#' plot(mod1, what = "IRF", IRF.args=list(item = 10,errorbar = TRUE))
#' # plot mastery probabilities for individuals 4 and 10
#' plot(mod1, what = "mp", mp.args = list(id = c(4,10)))
#'}

#' @export
plot.GDINA <-
  function(x, what = "IRF", IRF.args = list(item = "all", errorbar = FALSE, SE.type = 2),
           mp.args = list(id = 1, att.names = NULL), ...)
  {
    if(class(x)!="GDINA") stop("x must be of class GDINA.",call. = FALSE)

    if(toupper(what)=="IRF"){
      default.IRF.args = list(item = "all", errorbar = FALSE, SE.type = 2)
      IRF.args <- modifyList(default.IRF.args,IRF.args)
      if(extract(x,"sequential")){
        tit <- "Processing functions"
      }else{
        tit <- "Item success probabilities"
      }

      lc <- p <- upper <- lower <- NULL
      if (IRF.args$errorbar) se <- extract(x,what = "catprob.se",SE.type=IRF.args$SE.type)
      ip <- extract(x,what = "catprob.parm")
      if(IRF.args$item == "all"){
        item <- seq_len(length(ip))
      }else{
          item <- IRF.args$item
        }
      for (j in item){
        tmp.obj <- ip[[j]]

        tmp.name <- gsub("P\\(","",names(tmp.obj))
        tmp.name <- gsub("\\)","",tmp.name)
        # names(tmp.obj) <- tmp.name
        lc <- factor(tmp.name,levels = tmp.name)

        if(IRF.args$errorbar){
          lower=tmp.obj-se[[j]]
          lower[lower<0] <- 0
          upper=tmp.obj+se[[j]]
          upper[upper>1] <- 1
          dat <- data.frame(lc = lc,p = tmp.obj,lower=lower,upper=upper)
          print(ggplot2::ggplot(data = dat, aes(x = lc, y = p)) +
                  geom_bar(stat = "identity", position = "dodge") +
                  geom_errorbar(aes(ymax=upper,ymin=lower), position = "dodge", width = 0.15) +
                  ylim(0,1)+
                  labs(x = "Latent groups", y = "Probability of success",
                       title = paste(tit,"for", extract(x,"item.names")[j])))
        }else{
          dat <- data.frame(lc = lc,p = tmp.obj)
          print(ggplot2::ggplot(data = dat, aes(x = lc, y = p)) +
                  geom_bar(stat = "identity", position = "dodge") +
                  ylim(0,1)+
                  labs(x = "Latent groups", y = "Probability of success",
                       title = paste(tit,"for", extract(x,"item.names")[j])))
        }



      }
    }else if(tolower(what)=="mp"){
      default.mp.args = list(id = 1, type = "barchart", att.names = NULL)
      mp.args <- modifyList(default.mp.args,mp.args)
      df <- personparm(x,"mp")
      if(is.null(mp.args$att.names)){
        att.names <- colnames(df)
      }else{
        att.names <- mp.args$att.names
      }
      np <- length(mp.args$id)
      if(np>1){
        dff <- c(t(df[mp.args$id,]))
        dat <- data.frame(att = rep(att.names,np),mp = dff,id = factor(rep(mp.args$id,each = ncol(df))))
        print(ggplot2::ggplot(data = dat, ggplot2::aes_string(x = "att", y = "mp")) +
                geom_bar(stat = "identity", position = "dodge",ggplot2::aes_string(fill = "id")) +
                ylim(0,1)+
                labs(x = "Attributes", y = "Mastery probabilities",
                     title = paste("Mastery probabilities")))
      }else{
        dff <- c(df[mp.args$id,])
        dat <- data.frame(att = att.names,mp = dff,id = factor(rep(mp.args$id,ncol(df))))
        print(ggplot2::ggplot(data = dat, ggplot2::aes_string(x = "att", y = "mp")) +
                geom_bar(stat = "identity", position = "dodge") +
                ylim(0,1)+
                labs(x = "Attributes", y = "Mastery probabilities",
                     title = paste("Mastery probabilities for individual",mp.args$id)))
      }

    }


  }


#' Item fit plots
#'
#' Create plots of bivariate heatmap for item fit
#'
#' @param x model object of class \code{itemfit}
#' @param ... additional arguments
#' @seealso \code{\link{GDINA}}, \code{\link{itemfit}}
# #' @describeIn itemfit create bivariate heatmap plots
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#'
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' ift <- itemfit(fit)
#' plot(ift)
#'}
#' @export
plot.itemfit <- function(x,...){
  item.pair.1 <- item.pair.2 <- unadj.pvalue <- test.adj.pvalue <- NULL
  print(ggplot2::ggplot(extract.itemfit(x,"logOR"),
                        aes(x=factor(item.pair.2),
                            y=factor(item.pair.1),
                            fill=unadj.pvalue))+
          geom_tile()+ scale_fill_gradient(low="red",
                                           high="gray",
                                           limits=c(0,0.05))+
          theme_bw() +
          labs(x = "Items", y = "Items",
               title = "Heatmap plot for unadjusted p-values of log odds ratio"))

  print(ggplot2::ggplot(extract.itemfit(x,"logOR"),
                        aes(x=factor(item.pair.2),
                            y=factor(item.pair.1),
                            fill=test.adj.pvalue))+
          geom_tile()+ scale_fill_gradient(low="red",
                                           high="gray",
                                           limits=c(0,0.05))+
          theme_bw() +
          labs(x = "Items", y = "Items",
               title = "Heatmap plot for adjusted p-values of log odds ratio"))

  print(ggplot2::ggplot(extract.itemfit(x,"r"),
                        aes(x=factor(item.pair.2),
                            y=factor(item.pair.1),
                            fill=unadj.pvalue))+
          geom_tile()+ scale_fill_gradient(low="red",
                                           high="gray",
                                           limits=c(0,0.05))+
          theme_bw() +
          labs(x = "Items", y = "Items",
               title = "Heatmap plot for unadjusted p-values of transformed correlation"))

  print(ggplot2::ggplot(extract.itemfit(x,"r"),
                        aes(x=factor(item.pair.2),
                            y=factor(item.pair.1),
                            fill=test.adj.pvalue))+
          geom_tile()+ scale_fill_gradient(low="red",
                                           high="gray",
                                           limits=c(0,0.05))+
          theme_bw() +
          labs(x = "Items", y = "Items",
               title = "Heatmap plot for adjusted p-values of transformed correlation"))

}


#' Mesa plot for Q-matrix validation
#'
#' The mesa plot was first proposed by de la Torre and Ma (2016) for graphically illustrating the best q-vector(s) for each item.
#' The q-vector on the edge of the mesa is likely to be the best q-vector.
#'
#' @param x model object of class \code{Qvalidation}
#' @param item a vector specifying which item(s) the plots are drawn for
#' @param type types of the plot. It can be \code{"best"} or \code{"all"}. If \code{"best"},
#'     for all q-vectors requiring the same number of attributes, only the one with the largest PVAF
#'     is plotted, which means \eqn{K_j} q-vectors are plotted; If \code{"all"}, all q-vectors
#'     will be plotted.
#' @param eps the cutoff for PVAF. If not \code{NULL}, it must be a value between 0 and 1. A horizontal line will be drawn accordingly.
#' @param no.qvector the number of q vectors that need to be plotted when \code{type="all"}. The default is 10,
#'        which means the 10 q vectors with the largest PVAFs are plotted.
#' @param data.label logical; To show data label or not?
#' @param original.q.label logical; print the label showing the original q-vector or not?
#' @param auto.ylim logical; create y range automatically or not?
#' @param ... additional arguments passed to \code{plot} function
#' @seealso \code{\link{Qval}}, \code{\link{autoGDINA}}
#' @examples
#'\dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' Q[1,] <- c(0,1,0)
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' out <- Qval(mod1,eps = 0.9)
#' item <- c(1,2,10)
#' plot(out,item=item,data.label=FALSE,type="all")
#' plot(out,item=10,type="best",eps=0.95)
#' plot(out,item=10,type="all",no.qvector=5)
#'}
#'
#' @references
#'
#' de la Torre, J., & Ma, W. (2016, August). Cognitive diagnosis modeling: A general framework approach and its implementation in R. A Short Course at the Fourth Conference on Statistical Methods in Psychometrics, Columbia University, New York.
#'

#' @export
plot.Qval <-
  function(x, item, type = "best", no.qvector = 10,
           data.label = TRUE,eps = "auto",
           original.q.label = FALSE,auto.ylim = TRUE,...)
  {
    if(eps=="auto") eps <- x$eps
    Q <- extract.Qval(x,"Q")
    K <- ncol(Q)
    L <- (2^K-1) # L-1
    patt <- attributepattern(K)
    fullPVAF <- extract.Qval(x,"PVAF")
    if(tolower(type)=="all"){

      if (L<no.qvector) no.qvector <- L
      for (y in item){
        #which one is the true
        locy0 <- which(apply(patt[-1,],1,function(x){all(x==Q[y,])}))

        locy <- no.qvector-(L-which(order(fullPVAF[,y],decreasing = F)==locy0))
        ordered.PVAF.j <- sort(fullPVAF[,y],decreasing = FALSE)
        graphics::plot(ordered.PVAF.j[(L-no.qvector+1):L],xaxt="n",type="o",ylab = "PVAF",xlab="q-vectors",
                       main = paste("Mesa Plot for Item",y),ylim = c(0,1),...)
        axis(1,at=c(1:no.qvector),labels = names(ordered.PVAF.j[(L-no.qvector+1):L]))
        if (locy>0){
          points(locy,fullPVAF[locy0,y],col="red",pch=19)
        }
        if (!is.null(eps)&&eps>0&&eps<1) abline(h=eps,lty=3);text(1.5,eps+0.03,paste("eps =",eps))
        if (original.q.label) text(no.qvector-1,0.15,paste("original q-vector:\n",names(fullPVAF[,y])[locy0]))
        if (auto.ylim) ylim = c(max(0,round(min(ordered.PVAF.j)-0.1,1)),1) else ylim=c(0,1)
        yloc <- ordered.PVAF.j[(L-no.qvector+1):L]-diff(ylim)/15
        yloc[yloc<=ylim[1]] <- yloc[yloc<=ylim[1]] + 2 * diff(ylim)/15
        if (data.label) text(c(1:no.qvector),
                             yloc,
                             ordered.PVAF.j[(L-no.qvector+1):L])

      }
    }else if(tolower(type)=="best"){
      fullPVAF <- rbind(0,fullPVAF)
      Kj <- rowSums(patt)
      bestPVAF <- aggregate(fullPVAF,by=list(Kj),max)[,-1]
      # bestPVAF <- rbind(0,bestPVAF) # add 0s
      label.bestPVAF <- apply(patt,1,paste0,collapse = "")
      for(j in item){
        bestloc <- match(bestPVAF[,j],fullPVAF[,j])
        if (auto.ylim) ylim = c(max(0,round(min(bestPVAF[,j])-0.1,1)),1) else ylim=c(0,1)
        graphics::plot(bestPVAF[,j],xaxt="n",type="o",ylab = "PVAF",xlab="q-vectors",
                       main = paste("Mesa Plot for Item",j),ylim = ylim,...)
        graphics::axis(1,at=c(1:nrow(bestPVAF)),labels = label.bestPVAF[bestloc])
        if (!is.null(eps)&&eps>0&&eps<1) abline(h=eps,lty=3);text(1.5,eps+0.03,paste("eps =",eps))
        yloc <- bestPVAF[,j]-diff(ylim)/15
        yloc[yloc<=ylim[1]] <- yloc[yloc<=ylim[1]] + 2 * diff(ylim)/15
        if (data.label) graphics::text(c(1:nrow(bestPVAF)),yloc,bestPVAF[,j])
        locy0 <- which(apply(patt,1,function(x){
          all(x==Q[j,])}))
        if(locy0%in%bestloc) graphics::points(which(bestloc==locy0),fullPVAF[locy0,j],col="red",pch=19)
        if (original.q.label) text(K-1,ylim[1]+diff(ylim)/6,paste("original q-vector:\n",names(fullPVAF[,j])[locy0]))
      }
    }

  }


