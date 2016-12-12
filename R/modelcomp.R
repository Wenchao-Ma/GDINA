#' Item-level model comparison using Wald test
#'
#' This function evaluates whether the saturated G-DINA model can be replaced by reduced
#' CDMs without significant loss in model data fit for each item. See de la Torre and Lee (2013),
#' and Ma, Iaconangelo and de la Torre (2016) for details. This function also calculates the dissimilarity
#' between the reduced models and the G-DINA model, which can be viewed as a measure of effect size (Ma, Iaconangelo & de la Torre, 2016).
#'
#' @param GDINA.obj An estimated model object of class \code{GDINA}
#' @param item a vector of items to specify which items the Wald test is applied to
#' @param models a vector specifying which reduced CDMs are possible reduced CDMs for each
#'   item. The default is "DINA","DINO","ACDM","LLM",and "RRUM".
#' @param DS whether dissimilarity index should be calculated? \code{FALSE} is the default.
#' @param varcov Optional; user specified variance-covariance matrix. If supplied, it must
#'   be a list of length \eqn{J}, giving the variance covariance matrix of item success probability for each item.
#'   The default is \code{NULL}, in which case, the estimated variance-covariance matrix from the GDINA function
#'   is used.
#' @param SE.type the type of standard error estimates.
#'
#' @return an object of class \code{modelcomp}. Elements that can be
#' extracted using \code{extract} method include
#' \describe{
#' \item{wald}{wald statistics}
#' \item{wald.p}{p-values associated with the wald statistics}
#' \item{DS}{dissimilarity between G-DINA and other CDMs}
#' }
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#'
#' @seealso \code{\link{GDINA}}, \code{\link{autoGDINA}}
#'@references
#' de la Torre, J., & Lee, Y. S. (2013). Evaluating the wald test for item-level comparison of
#' saturated and reduced models in cognitive diagnosis. \emph{Journal of Educational Measurement, 50}, 355-373.
#'
#' Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity, model selection and attribute classification.
#' \emph{Applied Psychological Measurement, 40}, 200-217.
#'
#'@examples
#'\dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' # --- GDINA model ---#
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' mod1
#' wmod1 <- modelcomp(mod1)
#' wmod1
#' # wald statistics
#' extract(wmod1,"wald")
#' #p values
#' extract(wmod1,"wald.p")
#' wmod1r <- modelcomp(mod1,DS=TRUE)
#' #dissimilarity index
#' extract(wmod1r,"DS")
#' wmod2 <- modelcomp(mod1,models = c("DINA","DINO"))
#' wmod2
#' wmod3 <- modelcomp(mod1,item=c(8,9,10),DS=FALSE)
#' wmod3
#' }
#' @import MASS
#' @export
modelcomp <- function(GDINA.obj,item="all",DS=FALSE, SE.type = 2,
                           models=c("DINA","DINO","ACDM","LLM","RRUM"),
                           varcov = NULL){
if(internalextract(GDINA.obj,"sequential")) stop("Model selection is not available for sequential models.",call. = FALSE)
  if (any(internalextract(GDINA.obj,"models")!="GDINA")) stop ("Implementing the Wald test for model comparison requires all items are fitted by the G-DINA model.",call. = FALSE)
  models <- match.arg(models, several.ok = TRUE)
  Q <- internalextract(GDINA.obj,"Q")
  Kjs <- rowSums(Q)
  Ks <- cumsum(2^Kjs)
  W <- pvalues <- df <- NULL
  if(length(item)==1){
    if(tolower(item)=="all") {
    item <- which(rowSums(Q)>1)
  }else {
    if (!is.numeric(item)) stop("item must be specified as item numbers.",call. = FALSE)
    if (sum(Q[item,])==1) {
      item <- NULL
      stop("Wald test can only be used for the items requiring more than 1 attributes.",call. = FALSE)
    }
  }
    }else {
    if (!all(sapply(item,is.numeric))) stop("Items must be specified as item numbers.",call. = FALSE)

    item <- intersect(which(rowSums(Q)>1),item)
  }
  if (!is.null(item)){

    for (j in item){ # for each item
      wald <- p <- dfj <- 0
      Kj <- rowSums(Q[item,]) # Kj for each item
      RDA <- RDINA(unique(Kj))
      RDO <- RDINO(unique(Kj))
      RAM <- RACDM(unique(Kj))

      # variance covariance matrix for item j
      if (is.null(varcov)) {
        cov <- internalextract(GDINA.obj,"catprob.cov",type = SE.type)
        ind <- cov$index
        vcov <- cov$cov[ind[which(ind[,1]==j),2],ind[which(ind[,1]==j),2]]
        }
      else{vcov <- varcov[[j]]}

      if ("DINA" %in% models){
        w <- t(RDA[[Kjs[j]]]%*%internalextract(GDINA.obj,"catprob.parm")[[j]])%*%
          MASS::ginv(RDA[[Kjs[j]]]%*%vcov%*%t(RDA[[Kjs[j]]]))%*%
          (RDA[[Kjs[j]]]%*%internalextract(GDINA.obj,"catprob.parm")[[j]])
        wald <- c(wald,w)
        dfj <- c(dfj,2^Kjs[j]-2)
        p <- c(p,1-pchisq(w,2^Kjs[j]-2))
      }else{
        wald <- c(wald,NA)
        dfj <- c(dfj,NA)
        p <- c(p,NA)
      }
      if ("DINO" %in% models){
        w <- t(RDO[[Kjs[j]]]%*%internalextract(GDINA.obj,"catprob.parm")[[j]])%*%
          MASS::ginv(RDO[[Kjs[j]]]%*%vcov%*%t(RDO[[Kjs[j]]]))%*%
          (RDO[[Kjs[j]]]%*%internalextract(GDINA.obj,"catprob.parm")[[j]])
        wald <- c(wald,w)
        dfj <- c(dfj,2^Kjs[j]-2)
        p <- c(p,1-pchisq(w,2^Kjs[j]-2))
      }else{
        wald <- c(wald,NA)
        dfj <- c(dfj,NA)
        p <- c(p,NA)
      }
      if ("ACDM" %in% models){
        w <- t(RAM[[Kjs[j]]]%*%internalextract(GDINA.obj,"catprob.parm")[[j]])%*%
          MASS::ginv(RAM[[Kjs[j]]]%*%vcov%*%t(RAM[[Kjs[j]]]))%*%
          (RAM[[Kjs[j]]]%*%internalextract(GDINA.obj,"catprob.parm")[[j]])
        wald <- c(wald,w)
        dfj <- c(dfj,2^Kjs[j]-Kjs[j]-1)
        p <- c(p,1-pchisq(w,2^Kjs[j]-Kjs[j]-1))
      }else{
        wald <- c(wald,NA)
        dfj <- c(dfj,NA)
        p <- c(p,NA)
      }
      if ("LLM" %in% models){
        fpj <- logit(internalextract(GDINA.obj,"catprob.parm")[[j]])
        grad_fpj <- solve(diag(internalextract(GDINA.obj,"catprob.parm")[[j]]*(1-internalextract(GDINA.obj,"catprob.parm")[[j]])))
        var_fpj <- grad_fpj%*%vcov%*%t(grad_fpj)
        w <- t(RAM[[Kjs[j]]]%*%fpj)%*%
          MASS::ginv(RAM[[Kjs[j]]]%*%var_fpj%*%t(RAM[[Kjs[j]]]))%*%
          (RAM[[Kjs[j]]]%*%fpj)
        wald <- c(wald,w)
        dfj <- c(dfj,2^Kjs[j]-Kjs[j]-1)
        p <- c(p,1-pchisq(w,2^Kjs[j]-Kjs[j]-1))
      }else{
        wald <- c(wald,NA)
        dfj <- c(dfj,NA)
        p <- c(p,NA)
      }
      if ("RRUM" %in% models){
        fpj <- log(internalextract(GDINA.obj,"catprob.parm")[[j]])
        grad_fpj <- solve(diag(internalextract(GDINA.obj,"catprob.parm")[[j]]))
        var_fpj <- grad_fpj%*%vcov%*%t(grad_fpj)
        w <- t(RAM[[Kjs[j]]]%*%fpj)%*%
          MASS::ginv(RAM[[Kjs[j]]]%*%var_fpj%*%t(RAM[[Kjs[j]]]))%*%
          (RAM[[Kjs[j]]]%*%fpj)
        wald <- c(wald,w)
        dfj <- c(dfj,2^Kjs[j]-Kjs[j]-1)
        p <- c(p,1-pchisq(w,2^Kjs[j]-Kjs[j]-1))
      }else{
        wald <- c(wald,NA)
        dfj <- c(dfj,NA)
        p <- c(p,NA)
      }
      W <- rbind(W,wald)
      pvalues <- rbind(pvalues,p)
      df <- rbind(df,dfj)
    }
  }
  W <- W[,-1]
df <- df[,-1]
  pvalues <- pvalues[,-1]

  # W = round(W,digits)
  # pvalues = round(pvalues,digits)

  ##### dissimilarity index
ds.f <- NULL
if(DS){
  ds <- lapply(internalextract(GDINA.obj,"catprob.parm")[item],function(p){
    unlist(lapply(models,function(m) DS(p,m)$DS))
  })
  ds <- do.call(rbind,ds)
  ds.f <- matrix(NA,nrow(ds),5)
if(length(models)<5) {
  ds.f[,match(models,c("DINA","DINO","ACDM","LLM","RRUM"))] <- ds
}else{
  ds.f <- ds
}
  colnames(ds.f) <- c("DINA","DINO","ACDM","LLM","RRUM")
  # ds.f <- round(ds.f,digits)
}


  colnames(W) <- colnames(df) <- colnames(pvalues) <- c("DINA","DINO","ACDM","LLM","RRUM")
  rownames(W) <- rownames(df) <- rownames(pvalues) <- paste("Item",item)
  out <- list(wald=W,wald.p=pvalues,df=df,DS=ds.f,models = models)
  class(out) <- "modelcomp"
  return(out)

}
