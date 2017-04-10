#' @include GDINA.R
#' @title  Plot item success probability
#'
#' @description   Create plots of item/category success probability for each latent group
#'
#' @param object model object of class \code{\link{GDINA}} or \code{\link{dif}}
#' @param item a vector specifying which item(s) the plots are drawn for
#' @param errorbar add error bar to the plot?
#' @param ... additional arguments
#' @seealso \code{\link{GDINA}}, \code{\link{autoGDINA}}, \code{\link{dif}}
#' @export
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' #plot item response functions for item 10
#' plotIRF(mod1,10)
#' plotIRF(mod1,9, errorbar = TRUE)
#'}
plotIRF <- function(object, item, errorbar = FALSE, ...){
  UseMethod("plotIRF")
}
#' @export
plotIRF.GDINA <-
  function(object, item, errorbar = FALSE, SE.type = 2,...)
  {
    if(extract(object,"sequential")) stop("IRF plot is not available for sequential models.",call. = FALSE)
    lc <- p <- upper <- lower <- NULL
if(class(object)!="GDINA") stop("object must be of class GDINA.",call. = FALSE)
    if (errorbar) se <- extract(object,what = "catprob.se",SE.type=SE.type)
    for (j in item){
      tmp.obj <- extract(object,what = "catprob.parm")[[j]]

      tmp.name <- gsub("P\\(","",names(tmp.obj))
      tmp.name <- gsub("\\)","",tmp.name)
      # names(tmp.obj) <- tmp.name
      lc <- factor(tmp.name,levels = tmp.name)

      if(errorbar){
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
                     title = paste("Item success probabilities for Item", j)))
      }else{
        dat <- data.frame(lc = lc,p = tmp.obj)
        print(ggplot2::ggplot(data = dat, aes(x = lc, y = p)) +
                geom_bar(stat = "identity", position = "dodge") +
                ylim(0,1)+
                labs(x = "Latent groups", y = "Probability of success",
                     title = paste("Item success probabilities for Item", j)))
      }



    }

  }


#' @export
plotIRF.dif <-
  function(object, item, errorbar = FALSE, SE.type = 2, ...)
  {
    # if(class(object)!="dif") stop("object must be of class dif",call. = FALSE)
lc <- p <- DIFgroups <- NULL
if (errorbar) {
  se <- extract(object$mg.est,what = "catprob.se",SE.type=SE.type)
}
J <- 0.5*extract(object$mg.est,"ncat")

    for (j in item){
      if (j>J) stop("Item number must not beyond the range.",call. = FALSE)
      tmp.obj <- c(extract(object$mg.est,what = "catprob.parm")[[j]],
                    extract(object$mg.est,what = "catprob.parm")[[j+J]])
      tmp.name <- gsub("P\\(","",names(tmp.obj))
      tmp.name <- gsub("\\)","",tmp.name)

      lc <- factor(x=tmp.name,levels=unique(tmp.name))

      if (errorbar) {
        sej <- c(se[[j]],se[[j+J]])
        lower=tmp.obj-sej
        lower[lower<0] <- 0
        upper=tmp.obj+sej
        upper[upper>1] <- 1
        dat <- data.frame(lc = lc,p = tmp.obj,lower=lower,upper=upper,
                          DIFgroups = as.factor(rep(unique(object$group),each = 0.5*length(tmp.obj))))
        print(ggplot2::ggplot(dat,aes(lc,p,fill = DIFgroups))+
                geom_bar(stat = "identity", position=position_dodge())+
                geom_errorbar(aes(ymax=upper,ymin=lower), position=position_dodge(.9), width = 0.15) +
                ylim(0,1)+
                labs(x = "Latent groups", y = "Probability of success",
                     title = paste("Item success probabilities for Item",j)))
      }else{
        dat <- data.frame(lc = lc,p = tmp.obj,
                          DIFgroups = as.factor(rep(unique(object$group),each = 0.5*length(tmp.obj))))
        print(ggplot2::ggplot(dat,aes(lc,p,fill = DIFgroups))+
                geom_bar(stat = "identity", position = "dodge")+
                ylim(0,1)+
                labs(x = "Latent groups", y = "Probability of success",
                     title = paste("Item success probabilities for Item",j)))
      }

    }
  }

