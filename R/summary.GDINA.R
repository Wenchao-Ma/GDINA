#' @include GDINA.R autoGDINA.R modelcomp.R itemfit.R GDI.R dif.R
#' @export
#' @describeIn GDINA print summary information
summary.GDINA <-
  function(object, ...)
  {
    output <- list("Loglikelihood"=extract(object,"logLik"),
    "Deviance"=extract(object,"deviance"),
    "AIC"=extract(object,"AIC"),
    "AIC penalty"=2*extract(object,"npar"),
    "AIC penalty due to item parameters"=2*extract(object,"npar.item"),
    "AIC penalty due to population parameters"=2*extract(object,"npar.att"),
    "AICc"=extract(object,"AICc"),
    "AICc penalty"=2*extract(object,"npar")*(extract(object,"npar")+1)/(extract(object,"nobs")-extract(object,"npar")-1),
    "AICc penalty due to item parameters"=2*extract(object,"npar.item")*(extract(object,"npar.item")+1)/(extract(object,"nobs")-extract(object,"npar.item")-1),
    "AICc penalty due to population parameters"=2*extract(object,"npar.att")*(extract(object,"npar.att")+1)/(extract(object,"nobs")-extract(object,"npar.att")-1),
    "BIC"=extract(object,"BIC"),
    "BIC penalty"=log(extract(object,"nobs"))*extract(object,"npar"),
    "BIC penalty due to item parameters"=log(extract(object,"nobs"))*extract(object,"npar.item"),
    "BIC penalty due to population parameters"=log(extract(object,"nobs"))*extract(object,"npar.att"),
    "CAIC"=extract(object,"CAIC"),
    "CAIC penalty"=(log(extract(object,"nobs"))+1)*extract(object,"npar"),
    "CAIC penalty due to item parameters"=(log(extract(object,"nobs"))+1)*extract(object,"npar.item"),
    "CAIC penalty due to population parameters"=(log(extract(object,"nobs"))+1)*extract(object,"npar.att"),
    "SABIC"=extract(object,"SABIC"),
    "SABIC penalty"=log((extract(object,"nobs")+2)/24)*extract(object,"npar"),
    "SABIC penalty due to item parameters"=log((extract(object,"nobs")+2)/24)*extract(object,"npar.item"),
    "SABIC penalty due to population parameters"=log((extract(object,"nobs")+2)/24)*extract(object,"npar.att"),
    "Attribute Prevalence"=extract(object,"prevalence"),
    "ngroup"=extract(object,"ngroup"),
    "Number of parameters"=extract.GDINA(object,"npar"),
    "Number of estimated item parameters" = extract.GDINA(object,"npar.item"),
    "Number of fixed item parameters"=extract.GDINA(object,"npar.fixeditem"),
    "Number of population parameters"=extract.GDINA(object,"npar.att"))
    class(output) <- "summary.GDINA"
    output
  }


#' @export
#' @describeIn autoGDINA print summary information
summary.autoGDINA <-
  function(object, ...)
  {

    fit <- data.frame(npar=c(extract(object$GDINA1.obj,"npar"),
                             extract(object$GDINA2.obj,"npar"),
                             extract(object$CDM.obj,"npar")),
                      logLik=round(c(logLik(object$GDINA1.obj),
                                     logLik(object$GDINA2.obj),
                                     logLik(object$CDM.obj)),2),
                      deviance=c(deviance(object$GDINA1.obj),
                                 deviance(object$GDINA2.obj),
                                 deviance(object$CDM.obj)),
                      AIC=c(AIC(object$GDINA1.obj),
                            AIC(object$GDINA2.obj),
                            AIC(object$CDM.obj)),
                      BIC=c(BIC(object$GDINA1.obj),
                            BIC(object$GDINA2.obj),
                            BIC(object$CDM.obj)),
                      row.names = c("Initial GDINA","GDINA using validated Q-matrix","Final CDMs"))
output <- list(fit=fit,Qval=object$Qval.obj,finalmodel=object$CDM.obj)
    class(output) <- "summary.autoGDINA"
    output
  }

 #' @export
 #' @describeIn itemfit print summary information
 summary.itemfit <-
   function(object, ...)
   {
     cat("\nItem-level fit statistics\n")
     print(extract.itemfit(object,"maxitemfit"))
     invisible(extract.itemfit(object,"maxitemfit"))
   }

 #' @export
 #' @param object dif object for S3 method
 #' @describeIn dif print summary information
 summary.dif <-
   function(object, ...)
   {
     cat("\nDIF Analysis Results\n")
     cat("\nWald statistics and p-values:\n")
     print(object$test)
     if(!is.null(object$p.adjust.methods))
       cat("\np-value adjustment method:", object$p.adjust.methods, "\n")
     invisible(object$test)
   }

 #' @export
 #' @describeIn Qval print summary information
 summary.Qval <-
   function(object, ...)
   {
     cat("\nQ-matrix validation\n")
     cat("PVAF threshold - eps =", extract.Qval(object,"eps"))
     cat("\nUse extract to extract various elements.")
     invisible(NULL)

   }


 #' @export
 #' @describeIn modelcomp print summary information
 summary.modelcomp <-
   function(object, ...)
   {
     cat("\nItem-level model comparison:\n")
     cat("Test statistics for items requiring two or more attributes:\n")
     wald <- extract.modelcomp(object,"stats")
     print(wald[,colSums(is.na(wald))==0])
     cat("\np-values for items requiring two or more attributes:\n")
     p <- extract.modelcomp(object,"pvalues")
     print(p[,colSums(is.na(p))==0])
   }
