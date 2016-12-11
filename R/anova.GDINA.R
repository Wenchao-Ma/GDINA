#' @include GDINA.R
#' @title Model Comparison
#' @description To compare two \code{GDINA} models, use method \code{\link{anova}}.
#' @export
#' @aliases anova.GDINA
#' @note anova function does NOT check whether models compared are nested or not.
#' @describeIn GDINA Model comparison using likelihood ratio test
anova.GDINA <- function(object, ...)
  {
    obj <- match.call()
    objects <- list(object, ...)
    if (length(objects) != 2){
      stop("Only two models can be compared using anova.\n")
    }

    # print(objects)
    GDINA.obj.1 <- objects[[1]]
    GDINA.obj.2 <- objects[[2]]

    if(class(GDINA.obj.1)!="GDINA"|class(GDINA.obj.2)!="GDINA") stop("Inputs must be objects from class GDINA.",call. = FALSE)
    out <- data.frame(npar=c(internalextract(GDINA.obj.1,"npar"),internalextract(GDINA.obj.2,"npar")),
                      logLik=c(internalextract(GDINA.obj.1,"logLik"),internalextract(GDINA.obj.2,"logLik")),
                      deviance=c(internalextract(GDINA.obj.1,"deviance"),internalextract(GDINA.obj.2,"deviance")),
                      AIC=c(internalextract(GDINA.obj.1,"AIC"),internalextract(GDINA.obj.2,"AIC")),
                      BIC=c(internalextract(GDINA.obj.1,"BIC"),internalextract(GDINA.obj.2,"BIC")))
    # print(out)
    colnames(out) <- c("#par","logLik","deviance","AIC","BIC")
    rownames(out) <- paste(obj)[-1]
if(npar(GDINA.obj.1)[1]!=npar(GDINA.obj.2)[1]){
  delchi <- NULL
  pval <- round(pchisq(abs(internalextract(GDINA.obj.1,"deviance")-internalextract(GDINA.obj.2,"deviance")),
                       abs(npar(GDINA.obj.1)[1]-npar(GDINA.obj.2)[1]),lower.tail = FALSE),4)
  pval <- ifelse(pval<0.001,"<0.001",pval)
  LR <- ifelse(npar(GDINA.obj.1)[1]>npar(GDINA.obj.2)[1],
               internalextract(GDINA.obj.2,"deviance")-internalextract(GDINA.obj.1,"deviance"),
               internalextract(GDINA.obj.1,"deviance")-internalextract(GDINA.obj.2,"deviance"))
  if (LR<0) pval <- "NA"
  delchi <- data.frame(chisq=c("",round(LR,4)),
                       df=c("",round(abs(npar(GDINA.obj.1)[1]-npar(GDINA.obj.2)[1]),3)),
                       pvalue=c("",pval))
  colnames(delchi) <- c("chisq","df","p-value")
  out <- cbind(out,delchi)

}
    return(out)
  }
