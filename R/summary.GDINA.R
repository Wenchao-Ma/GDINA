#' @include GDINA.R autoGDINA.R modelcomp.R itemfit.R GDI.R dif.R
#' @export
#' @describeIn GDINA print summary information
summary.GDINA <-
  function(object, ...)
  {
    cat("\nTest Fit Statistics\n\n")
    cat("Loglikelihood =", format(round(extract.GDINA(object,"logLik"), 2), nsmall = 2), "\n")
    cat("Deviance      =", format(round(extract.GDINA(object,"deviance"), 2), nsmall = 2), "\n")
    cat("AIC           =", format(round(extract.GDINA(object,"AIC"), 2), nsmall = 2),"\n")
    cat("AIC Penalty   =",2*extract.GDINA(object,"npar"),"\n")
    cat("  AIC penalty due to item parameters        =", 2*extract.GDINA(object,"npar.item"), "\n")
    cat("  AIC penalty due to population parameters  =", 2*extract.GDINA(object,"npar.att"), "\n")
    cat("BIC           =", format(round(extract.GDINA(object,"BIC"), 2), nsmall = 2),"\n")
    cat("BIC penalty   =",format(round(log(extract.GDINA(object,"nobs"))*extract.GDINA(object,"npar"),2), nsmall = 2),"\n")
    cat("  BIC penalty due to item parameters        =", format(round(log(extract.GDINA(object,"nobs"))*extract.GDINA(object,"npar.item"),2), nsmall = 2), "\n")
    cat("  BIC penalty due to population parameters  =", format(round(log(extract.GDINA(object,"nobs"))*extract.GDINA(object,"npar.att"),2), nsmall = 2), "\n")
    if(extract.GDINA(object,"ngroup")==1){
      cat("\nAttribute Prevalence\n\n")
      print(round(extract.GDINA(object,"prevalence"),extract.GDINA(object,"digits")))
    }

    if( extract.GDINA(object,"natt")<4){
    cat("\nPosterior Weights\n\n")
    print(extract.GDINA(object,"posterior.prob"))
    }
    invisible(NULL)
  }


#' @export
#' @describeIn autoGDINA print summary information
summary.autoGDINA <-
  function(object, ...)
  {
    cat("\n")
    out <- data.frame(npar=c(extract(object$GDINA1.obj,"npar"),
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

    print(out)
    cat("\n")
    if(object$options$Qvalid){
      print(object$Qval.obj)
    }
    if(object$options$modelselection){
      cat("\nSelected Models:\n\n")
      print(extract.GDINA(object$CDM.obj,"models"))
    }
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
 #' @describeIn dif print summary information
 summary.dif <-
   function(object, ...)
   {
     cat("\nItem success probabilities for two groups\n")
     out <- mapply(rbind,extract.GDINA(object$CDM1,what = "itemprob.parm"),
                   extract.GDINA(object$CDM2,what = "itemprob.parm"),SIMPLIFY = F)
     out <- lapply(out,function(x){rownames(x) <- c("Group1.Est.","Group2.Est.");x})
     print(out)
     invisible(out)
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
     cat("Wald statistics for items requiring two or more attributes:\n")
     wald <- extract.modelcomp(object,"wald")
     print(wald[,colSums(is.na(wald))==0])
     cat("\nWald test p-values for items requiring two or more attributes:\n")
     p <- extract.modelcomp(object,"wald.p")
     print(p[,colSums(is.na(p))==0])
   }
