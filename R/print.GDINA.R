#' @include GDINA.R autoGDINA.R modelcomp.R itemfit.R GDI.R dif.R s3GDINA.R
#' @export
print.GDINA <-
  function(x, ...)
  {
    cat("Call:\n", paste(deparse(extract.GDINA(x,"call")), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    packageinfo <- utils::packageDescription("GDINA")
    cat( paste( "  GDINA version " , packageinfo$Version , " (" , packageinfo$Date , ")" , sep="") , "\n" )
    cat("===============================================\n")
    cat("Data\n")
    cat("-----------------------------------------------\n")
    cat("# of individuals    groups    items         \n")
    cat("    ",sprintf("%11d",extract.GDINA(x,"nobs")),sprintf("%9d",extract.GDINA(x,"ngroup")),sprintf("%8d",extract.GDINA(x,"nitem")))
    # cat("Number of items       =", extract.GDINA(x,"nitem"), "\n")
    # cat("Number of individuals =", extract.GDINA(x,"nobs"), "\n")
    # cat("Number of attributes  =", extract.GDINA(x,"natt"), "\n")
    # cat("Number of groups      =", extract.GDINA(x,"ngroup"), "\n")
    M <- c("User-defined Function","GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM")
    # cat("Response level        =",ifelse(max(extract.GDINA(x,"dat"),na.rm = TRUE)>1,"Polytomous","Dichotomous"),"\n")
    cat("\n===============================================\n")
    # cat("\n-----------------------------------------------\n")
    cat("Model")
    cat("\n-----------------------------------------------\n")
    lv <- NULL
    if(extract(x,"latent.var")=="bugs"){
      lv <- "Bug"
    }
    cat("Fitted model(s)       =", lv,unique(extract.GDINA(x,"models")), "\n")
    cat("Attribute structure   =",extract(x,"att.dist"),"\n")
    if (any(extract.GDINA(x,"att.dist")=="higher.order")) cat("Higher-order model    =",extract(x,"higher.order")$model,"\n")
    tmp <- max(extract.GDINA(x,"Q"))
    cat("Attribute level       =",ifelse(tmp>1,"Polytomous","Dichotomous"),"\n")
    cat("===============================================\n")
    cat("Estimation\n")
    cat("-----------------------------------------------\n")
    cat("Number of iterations  =", max(extract.GDINA(x,"nitr")), "\n")
    cat("For the final iteration:\n")
    cat("  Max abs change in item success prob. =", formatC(extract(x,"dif.p"), digits = 4, format = "f"), "\n")
    cat("  Max abs change in mixing proportions =", formatC(extract(x,"dif.prior"), digits = 4, format = "f"), "\n")
    cat("  Change in -2 log-likelihood          =", formatC(extract(x,"dif.LL"), digits = 4, format = "f"), "\n")
    cat("Time used             =", format(extract(x,"time"), digits = 4), "\n")
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
print.CA <-
  function(x, ...)
  {
    cat("Classification Accuracy \n")
    cat("\nTest level accuracy = ", round(x$tau,4), "\n")
    cat("\nPattern level accuracy: \n\n")
    print(round(x$tau_l,4))
    cat("\nAttribute level accuracy: \n\n")
    print(round(x$tau_k,4))
  }
#' @export
print.modelcomp <- function(x, ...)
{
  cat("\n",x$method,"statistics for items requiring two or more attributes:\n")
  ret <- extract.modelcomp(x,"stats")
  print(ret[,colSums(is.na(ret))==0],drop=FALSE)
  cat("\np-values for items requiring two or more attributes:\n")
  p <- extract.modelcomp(x,"pvalues")
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
print.modelfit <-
  function(x, ...)
  {
    cat("Model Fit Evaluation\n\n")
    cat("-2 log likelihood = ",round(x$Deviance,4)," number of parameters = ", x$npar,"\n")
    cat("AIC = ", round(x$AIC,4)," BIC = ", round(x$BIC,4),"\n")
    if(!is.null(x$M2)){
      cat("M2 = ", round(x$M2,4)," df = ", x$M2.df," p = ", round(x$M2.pvalue,4),"\n")
      cat("RMSEA = ", round(x$RMSEA,4)," with ",x$CI*100,"% CI: [",round(x$RMSEA.CI[1],4),",",round(x$RMSEA.CI[2],4),"]\n")

    }
    cat("SRMSR = ", round(x$SRMSR,4))
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
  cat("Loglik =", formatC(x$Loglikelihood,digits = 2, format = "f"), "\n")
  # cat("Deviance      =", formatC(x$Deviance,digits = 2, format = "f"), "\n")
  cat("AIC    =", formatC(x$AIC,digits = 2, format = "f")," | penalty   =",x$`AIC Penalty`,"\n")
  cat("BIC    =", formatC(x$BIC,digits = 2, format = "f")," | penalty   =",formatC(x$`BIC penalty`,digits = 2, format = "f"),"\n")
  cat("# par  =",formatC(x$`Number of parameters`,digits = 0, format = "d"), "\n")
    cat("\nAttribute Prevalence\n\n")
    ap <- lapply(x$`Attribute Prevalence`,round,digits=4)
    if(x$ngroup==1) {
      print(ap[[1]])
    }else{
      print(ap)
    }


  # cat("\nPosterior Weights\n\n")
  #
  #     print(round(x$`Posterior Weights`,4))

}

#'@export
print.anova.GDINA <- function(x,...){
  cat("\nInformation Criteria and Likelihood Ratio Test\n\n")
   LR <- x$LR
   for(m in 1:nrow(LR)) {
     if(LR$pvalue[m]<0.001) LR$pvalue[m] <- "<0.001"
     if(LR$df[m]==0) LR$chisq[m] <- LR$pvalue[m] <- LR$df[m] <- ""
     if(LR$chisq[m]<0) LR$df[m] <- LR$pvalue[m] <- ""
   }
   out <- cbind(x$IC,LR)
   colnames(out) <- c("#par","logLik","Deviance","AIC","BIC","chisq","df","p-value")
   print(out)
  if(nrow(x$IC)>2) cat("\nNotes: In LR tests, models were tested against",rownames(x$IC)[which.max(x$IC$npar)],"\n       LR test(s) do NOT check whether models are nested or not.")
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
print.npar.GDINA <- function(x,...){
  cat("No. of total parameters =",x$`No. of parameters`,"\n")
  cat("No. of item parameters =",x$`No. of item parameters`,"\n")
  cat("No. of population parameters =",x$`No. of population parameters`,"\n")
}
