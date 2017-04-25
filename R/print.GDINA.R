#' @include GDINA.R autoGDINA.R modelcomp.R itemfit.R GDI.R dif.R s3GDINA.R
#' @export
print.GDINA <-
  function(x, ...)
  {
    cat("\nThe Generalized DINA Model Framework  \n")
    packageinfo <- utils::packageDescription("GDINA")
    cat( paste( "   Version " , packageinfo$Version , " (" , packageinfo$Date , ")" , sep="") , "\n" )
    cat("\nCall:\n", paste(deparse(extract.GDINA(x,"call")), sep = "\n", collapse = "\n"),
        "\n", sep = "")
    cat("\nNumber of items       =", extract.GDINA(x,"nitem"), "\n")
    cat("Number of individuals =", extract.GDINA(x,"nobs"), "\n")
    cat("Number of attributes  =", extract.GDINA(x,"natt"), "\n")
    cat("Number of groups      =", extract.GDINA(x,"ngroup"), "\n")
    M <- c("GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM")
    cat("Number of iterations  =", extract.GDINA(x,"nitr"), "\n")
    cat("Fitted model(s)       =", unique(extract.GDINA(x,"models")), "\n")
    cat("Attribute structure   =",extract(x,"att.dist"),"\n")
    if (extract.GDINA(x,"ngroup")==1&&extract.GDINA(x,"att.dist")=="higher.order") cat("Higher-order model    =",extract.GDINA(x,"higher.order.model"),"\n")
    tmp <- max(extract.GDINA(x,"Q"))
    cat("Attribute level       =",ifelse(tmp>1,"Polytomous","Dichotomous"),"\n")
    cat("Response level        =",ifelse(max(extract.GDINA(x,"dat"),na.rm = TRUE)>1,"Polytomous","Dichotomous"),"\n")
    cat("\nNumber of parameters  =", extract.GDINA(x,"npar"), "\n")
    cat("  No. of item parameters       =",extract.GDINA(x,"npar.item"),"\n")
    cat("  No. of population parameters =",extract.GDINA(x,"npar.att"),"\n")
    cat("\nFor the last iteration:\n")
    cat("  Max abs change in success prob. =", formatC(extract(x,"dif.p"), digits = 4, format = "f"), "\n")
    cat("  Abs change in deviance          =", formatC(extract(x,"dif.LL"), digits = 4, format = "f"), "\n")
    cat("\nTime used             =", format(extract(x,"time"), digits = 4), "\n")
  }
#' @export
print.simGDINA <-
  function(x, ...)
  {
    cat("Data simulation using GDINA package \n")
    cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")
    cat("\nNumber of items       =", ncol(x$dat), "\n")
    cat("Number of individuals =", nrow(x$dat), "\n")
    cat("To extract components, use the method extract. \n")

  }
#' @export
print.modelcomp <- function(x, ...)
{
  cat("Wald statistics for items requiring two or more attributes:\n")
  wald <- extract.modelcomp(x,"wald")
  print(wald[,colSums(is.na(wald))==0],drop=FALSE)
  cat("\nWald test p-values for items requiring two or more attributes:\n")
  p <- extract.modelcomp(x,"wald.p")
  print(p[,colSums(is.na(p))==0])
}

#' @export
print.itemfit <-
  function(x, ...)
  {
    cat("Summary of Item Fit Analysis\n")
    cat("\nCall:\n", paste(deparse(x$options$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")
    cat("\n")
    p <- extract.itemfit(x,"p")
    r <- extract.itemfit(x,"r")
    logOR <- extract.itemfit(x,"logOR")
    testlevel.itemfit <- data.frame(p=c(mean(p$pstat[is.finite(p$pstat)],na.rm = TRUE),
                                        max(p$pstat[is.finite(p$pstat)],na.rm = TRUE),
                                        max(p$zstat[is.finite(p$zstat)],na.rm = TRUE),
                                        p$unadj.pvalue[which(p$zstat==max(p$zstat[is.finite(p$zstat)],na.rm = TRUE))],
                                        p$test.adj.pvalue[which(p$zstat==max(p$zstat[is.finite(p$zstat)],na.rm = TRUE))]),
                                    r=c(mean(r$rstat[is.finite(r$rstat)],na.rm = TRUE),
                                        max(r$rstat[is.finite(r$rstat)],na.rm = TRUE),
                                        max(r$zstat[is.finite(r$zstat)],na.rm = TRUE),
                                        r$unadj.pvalue[which(r$zstat==max(r$zstat[is.finite(r$zstat)],na.rm = TRUE))],
                                        r$test.adj.pvalue[which(r$zstat==max(r$zstat[is.finite(r$zstat)],na.rm = TRUE))]),
                                    l=c(mean(logOR$lstat[is.finite(logOR$lstat)],na.rm = TRUE),
                                        max(logOR$lstat[is.finite(logOR$lstat)],na.rm = TRUE),
                                        max(logOR$zstat[is.finite(logOR$zstat)],na.rm = TRUE),
                                        logOR$unadj.pvalue[which(logOR$zstat==max(logOR$zstat[is.finite(logOR$zstat)],na.rm = TRUE))],
                                        logOR$test.adj.pvalue[which(logOR$zstat==max(logOR$zstat[is.finite(logOR$zstat)],na.rm = TRUE))]))
    colnames(testlevel.itemfit) <- c("Proportion correct","Transformed correlation","Log odds ratio")
    rownames(testlevel.itemfit) <- c("mean[stats]","max[stats]",
                                     "max[z.stats]","p-value","adj.p-value")
    print(t(round(testlevel.itemfit,extract.itemfit(x,"digits"))))
    cat("Note: p-value and adj.p-value are associated with max[z.stats].")
    cat("\n      adj.p-values are based on the", extract.itemfit(x,"p.adjust.method"),"method.")
    if(any(is.na(p))|any(is.infinite(unlist(p)))) warning("Proportions have NA or Inf - check results!",call. = FALSE)
    if(any(is.na(r))|any(is.infinite(unlist(r)))) warning("Transformed correlations have NA or Inf - check results!",call. = FALSE)
    if(any(is.na(logOR))|any(is.infinite(unlist(logOR)))) warning("Log odds ratios have NA or Inf - check results!",call. = FALSE)
  }

#' @export
print.Qval <-
  function(x, ...)
  {
    sugQ <- data.frame(extract.Qval(x,"sug.Q"))
    oriQ <- data.frame(extract.Qval(x,"Q"))
    if(any(sugQ!=oriQ)){
      cat("Suggested Q-matrix: \n")
      sugQ[sugQ!=oriQ] <- paste0(sugQ[sugQ!=oriQ],"*")
      print(sugQ,right = FALSE)
      cat("Note: * denotes a modified element.\n")
    }else{
      cat("\nNo Q-matrix modifications are suggested.\n")
    }

  }

#' @export
print.dif <-
  function(x, ...)
  {
    cat("\nDifferential Item Functioning Detection\n")
    print(x$test)
cat("\nNote: adjusted pvalues are based on the",x$p.adjust.methods,"correction.\n")
  }

#' @export
print.autoGDINA <-
  function(x, ...)
  {
    packageinfo <- utils::packageDescription("GDINA")
    cat( paste( "\nGDINA Version " , packageinfo$Version , " (" , packageinfo$Date , ")" , sep="") , "\n" )
    cat("\nCall:\n", paste(deparse(x$options$ASEcall), sep = "\n", collapse = "\n"),
        "\n", sep = "")
    cat("\nInitial GDINA calibration [GDINA1.obj] -> ")
    if(x$options$Qvalid) {
      cat("\nQ-matrix validation [Qval.obj] -> ")
      cat("\nSecond GDINA calibration [GDINA2.obj] -> ")
    }else{
      cat("\nNo Q-matrix validation [Qval.obj = NULL; GDINA2.obj = GDINA1.obj] -> ")
    }
    if(x$options$modelselection) {
      cat("\nItem-level model selection [Wald.obj] -> \nSelected CDMs calibration [CDM.obj]\n\n")
    }else{
      cat("\nNo Item-level model selection [Wald.obj = NULL] -> \nFinal calibration [CDM.obj]\n\n")
    }

    if(x$options$Qvalid){
      cat(" - Q-matrix validation is based on eps = ",x$options$eps,"\n\n")
    }else{
      #cat("\n - Q-matrix validation is disabled\n\n")
    }
    if(x$options$modelselection){
      cat(" - Model selection is based on the selection rule of",x$options$modelselectionrule,"\n")
      cat(" - Reduced models include",x$options$reducedCDM,"\n")
    cat(" - The alpha level for the Wald test is",x$options$alpha.level,"\n")
    }else{
      #cat(" - Model selection is disabled\n")
    }



  }

#'@export
print.summary.GDINA <- function(x,...){
  cat("\nTest Fit Statistics\n\n")
  cat("Loglikelihood =", formatC(x$Loglikelihood,digits = 4, format = "f"), "\n")
  cat("Deviance      =", formatC(x$Deviance,digits = 4, format = "f"), "\n")
  cat("AIC           =", formatC(x$AIC,digits = 4, format = "f"),"\n")
  cat("AIC Penalty   =",x$`AIC Penalty`,"\n")
  cat("  AIC penalty due to item parameters        =", x$`AIC penalty due to item parameters`, "\n")
  cat("  AIC penalty due to population parameters  =", x$`AIC penalty due to population parameters`, "\n")
  cat("BIC           =", formatC(x$BIC,digits = 4, format = "f"),"\n")
  cat("BIC penalty   =",formatC(x$`BIC penalty`,digits = 4, format = "f"),"\n")
  cat("  BIC penalty due to item parameters        =", formatC(x$`BIC penalty due to item parameters`,digits = 4, format = "f"), "\n")
  cat("  BIC penalty due to population parameters  =", formatC(x$`BIC penalty due to population parameters`,digits = 4, format = "f"), "\n")

    cat("\nAttribute Prevalence\n\n")
    ap <- lapply(x$`Attribute Prevalence`,round,digits=4)
    if(x$ngroup==1) {
      print(ap[[1]])
    }else{
      print(ap)
    }


  cat("\nPosterior Weights\n\n")

      print(round(x$`Posterior Weights`,4))

}

#'@export
print.anova.GDINA <- function(x,...){
  cat("\nInformation Criteria and Likelihood Ratio Tests\n\n")
   LR <- x$LR
   for(m in 1:nrow(LR)) {
     if(LR$pvalue[m]<0.001) LR$pvalue[m] <- "<0.001"
     if(LR$df[m]==0) LR$chisq[m] <- LR$pvalue[m] <- LR$df[m] <- ""
     if(LR$chisq[m]<0) LR$df[m] <- LR$pvalue[m] <- ""
   }
   out <- cbind(x$IC,LR)
   colnames(out) <- c("#par","logLik","Deviance","AIC","BIC","chisq","df","p-value")
   print(out)
  if(nrow(x$IC)>2) cat("\nNotes: In LR tests, models were tested against",rownames(x$IC)[which.max(x$IC$npar)],"\n       LR tests did NOT check whether models are nested or not.")
}

#'@export
print.summary.autoGDINA <- function(x,...){
  cat("\nAutomatic GDINA Analysis\n\n")
  cat("Model Data fit\n")
  print(x$fit)

  if(!is.null(x$Qval)){
    cat("\n")
    # cat("\nSuggested Q-matrix:\n\n")
    print(x$Qval)
  }
  if(!is.null(x$finalmodel)){
    cat("\nSelected Models:\n\n")
    print(extract(x$finalmodel,"models"))
  }

}


#'@export
print.AIC.GDINA <- function(x,...){
  cat("AIC =",sprintf("%.4f", x))
}
#'@export
print.BIC.GDINA <- function(x,...){
  cat("BIC =",sprintf("%.4f", x))
}
#'@export
print.logLik.GDINA <- function(x,...){
  cat("logLik =",sprintf("%.4f", x))
}
#'@export
print.deviance.GDINA <- function(x,...){
  cat("deviance =",sprintf("%.4f", x))
}
#'@export
print.npar.GDINA <- function(x,...){
  cat("No. of total parameters =",x$`No. of parameters`,"\n")
  cat("No. of item parameters =",x$`No. of item parameters`,"\n")
  cat("No. of population parameters =",x$`No. of population parameters`,"\n")
}
