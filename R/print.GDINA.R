#' @include GDINA.R autoGDINA.R modelcomp.R itemfit.R GDI.R dif.R
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
    if(extract.GDINA(x,"att.str")){
      strc <- "User specified"
    }else {
      if(extract.GDINA(x,"higher.order")){
        strc <- "Higher-order"
      }else{
        strc <- "Saturated"
      }
    }

    cat("Attribute structure   =",strc,"\n")
    if (extract.GDINA(x,"higher.order")) cat("Higher-order model    =",extract.GDINA(x,"higher.order.model"),"\n")
    tmp <- ifelse(extract.GDINA(x,"sequential"),max(extract.GDINA(x,"Q")),max(extract.GDINA(x,"Q")[,-c(1:2)]))
    cat("Attribute level       =",ifelse(tmp>1,"Polytomous","Dichotomous"),"\n")
    cat("Response level        =",ifelse(max(extract.GDINA(x,"dat"),na.rm = TRUE)>1,"Polytomous","Dichotomous"),"\n")
    cat("\nNumber of parameters  =", extract.GDINA(x,"npar"), "\n")
    cat("  No. of item parameters       =",extract.GDINA(x,"npar.item"),"\n")
    cat("  No. of population parameters =",extract.GDINA(x,"npar.att"),"\n")
    cat("\nFor the last iteration:\n")
    cat("  Max abs change in success prob. =", format(round(extract.GDINA(x,"dif.p"), 5),scientific = FALSE), "\n")
    cat("  Abs change in deviance          =", format(round(extract.GDINA(x,"dif.LL"), 2),scientific = FALSE), "\n")
    cat("\nTime used             =", format(round(extract.GDINA(x,"time"), 4),scientific = FALSE), "\n")

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
  print(wald[,colSums(is.na(wald))==0])
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

  }

#' @export
print.autoGDINA <-
  function(x, ...)
  {
    packageinfo <- utils::packageDescription("GDINA")
    cat( paste( "\nGDINA Beta Version " , packageinfo$Version , " (" , packageinfo$Date , ")" , sep="") , "\n" )
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


