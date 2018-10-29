library(GDINA)
shinyServer(function(input, output) {

  ######## INPUTS

  output$contents1 <- shiny::renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    x <- read.csv(inFile$datapath, header = input$header,
             sep = input$sep, quote = input$quote)
    head(x)
  })

  output$contents2 <- shiny::renderTable({
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
    y <- read.csv(inFile$datapath, header = input$header2,
             sep = input$sep2, quote = input$quote)
    head(y)
  })


  ##################
  #  Model Estimation
  ##################

  est.result <- eventReactive(input$goButton, {
    withProgress(message = 'Model Estimating', value = 0.9, {
    inFile1 <- input$file1
    dat <- read.csv(inFile1$datapath, header = input$header,
             sep = input$sep, quote = input$quote)

    inFile2 <- input$file2
    Q <- read.csv(inFile2$datapath, header = input$header2,
                    sep = input$sep2, quote = input$quote)
    hom <- NULL
    if(input$attdis==0){
      HOdist <- "saturated"
    }else if(input$attdis==1){
      HOdist <- "higher.order"
      hom <- "Rasch"
    }else if(input$attdis==2){
      HOdist <- "higher.order"
      hom <- "1PL"
    }else if(input$attdis==3){
      HOdist <- "higher.order"
      hom <- "2PL"
    }else if(input$attdis==4){
      HOdist <- "fixed"
    }

    if(input$type!="UM"){
      m <- input$type
    }else{
      m <- strsplit(input$mv,",")[[1]]
      # inFile3 <- input$usermodels
      # m <- scan(inFile3$datapath,character(),sep = ",")
    }

      est <- GDINA::GDINA(dat = dat, Q = Q, model = m,
                          verbose = 0,att.dist = HOdist,
                          higher.order = list(model = hom),
                          sequential = input$seq,
                          mono.constraint = input$mono)


    est
  })})


  #### generate menu

  output$summary <- shinydashboard::renderMenu({
    if (input$goButton == 0)
      return()
    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("Estimation Summary", icon = icon("info"), tabName = "summary")
    )
  })
  output$fit <- shinydashboard::renderMenu({
    if (input$goButton == 0)
      return()
    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("Model Fit", icon = icon("check-square-o"), tabName = "fit")
    )
  })
  output$par <- shinydashboard::renderMenu({
    if (input$goButton == 0)
      return()
    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("Parameter Estimates", icon = icon("superscript"), tabName = "par")
    )
  })
  output$qv <- shinydashboard::renderMenu({
    if (input$goButton == 0||input$qvalcheck == 0)
      return()

    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("Q-matrix Validation Outputs", icon = icon("th"), tabName = "Qval")
    )
  })
  output$msec <- shinydashboard::renderMenu({
    if (input$goButton == 0||input$modelsel == 0)
      return()
    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("Model selection Outputs", icon = icon("list"), tabName = "ms")
    )
  })

  output$menuplot <- shinydashboard::renderMenu({
    if (input$goButton == 0)
      return()
    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("Plots", icon = icon("bar-chart"), tabName = "plot")
    )
  })





  ##################
  # Summary
  ##################
  info <- shiny::reactive({
    summary.info <- function(object){
      cat("\nLoglikelihood =",extract(object,"logLik"))
      cat("\nDeviance      =",extract(object,"deviance"))
      cat("\nAIC           =",extract(object,"AIC"))
      cat("\n  AIC Penalty =",2*extract(object,"npar"))
      # cat("\n  AIC penalty due to item parameters=",2*extract(object,"npar.item"))
      # cat("\n  AIC penalty due to population parameters=",2*extract(object,"npar.att"))
      cat("\nBIC           =",extract(object,"BIC"))
      cat("\n  BIC penalty =",round(log(extract(object,"nobs"))*extract(object,"npar"),2))
      # cat("\n  BIC penalty due to item parameters=",log(extract(object,"nobs"))*extract(object,"npar.item"))
      # cat("\n  BIC penalty due to population parameters=",log(extract(object,"nobs"))*extract(object,"npar.att"))
    }
    summary.info(est.result())
  })
  iter.info <- shiny::reactive({
    est.info <- function(x) {
      # cat("\nThe Generalized DINA Model Framework  \n")
      # packageinfo <- utils::packageDescription("GDINA")
      # cat( paste( "   GDINA Version " , packageinfo$Version , " (" , packageinfo$Date , ")" , sep="") , "\n" )
      # cat(  "   Wenchao Ma & Jimmy de la Torre \n" )
      # cat("See https://wenchao-ma.github.io/GDINA for more information.\n")

      cat("\nNumber of items       =", extract(x,"nitem"), "\n")
      cat("Number of individuals =", extract(x,"nobs"), "\n")
      cat("Number of attributes  =", extract(x,"natt"), "\n")
      cat("Number of iterations  =", extract(x,"nitr"), "\n")
      cat("Fitted models         =", unique(c(extract(x,"models"))),"\n")


      tmp <- ifelse(extract(x,"sequential"),max(extract(x,"Q")),max(extract(x,"Q")[,-c(1:2)]))
      cat("Attribute level       =",ifelse(tmp>1,"Polytomous","Dichotomous"),"\n")
      cat("Response level        =",ifelse(max(extract(x,"dat"),na.rm = TRUE)>1,"Polytomous","Dichotomous"),"\n")
      cat("\nNumber of parameters  =", extract(x,"npar"), "\n")
      cat("  No. of item parameters       =",extract(x,"npar.item"),"\n")
      cat("  No. of population parameters =",extract(x,"npar.att"),"\n")
      cat("\nFor the last iteration:\n")
      cat("  Max abs change in success prob. =", format(round(extract(x,"dif.p"), 5),scientific = FALSE), "\n")
      cat("  Abs change in deviance          =", format(round(extract(x,"dif.LL"), 2),scientific = FALSE), "\n")
      cat("\nTime used             =", format(round(extract(x,"time"), 4),scientific = FALSE), "\n")

    }
    est.info(est.result())
  })

  iter.info2 <- shiny::reactive({
    est.info2 <- function(x) {
      # cat("\nAttribute prevalence\n")
      # print(round(extract(x,"prevalence")$`1`,4))
      # cat("\n\n")
      print(CA(x))
    }
    est.info2(est.result())
  })


  output$info <- shiny::renderPrint({
    info()
  })

  output$iter.info <- shiny::renderPrint({
    iter.info()
  })
  output$iter.info2 <- shiny::renderPrint({
    iter.info2()
  })

  itf <- shiny::reactive({
    fitcheck <- function(object){
      x <- itemfit(object)
      if(all(extract(object,"models_numeric")>=0)&&all(extract(object,"models_numeric")<=5)&&input$attdis==0){
        z <- modelfit(object)

        cat("\nM2=",z$M2,"( df=",z$M2.df,")","p-value=",round(z$M2.pvalue,4))
        cat("\nRMSEA = ", round(z$RMSEA,4)," with ",z$CI*100,"% CI: [",round(z$RMSEA.CI[1],4),",",round(z$RMSEA.CI[2],4),"]")
        cat("\nSRMSR = ", round(z$SRMSR,4),"\n\n")
      }

      p <- extract(x,"p")
      r <- extract(x,"r")
      logOR <- extract(x,"logOR")
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
      print(t(round(testlevel.itemfit,4)))
      cat("Note: p-value and adj.p-value are associated with max[z.stats].")
      cat("\n      adj.p-values are based on the", extract(x,"p.adjust.method"),"method.")
      if(any(is.na(p))|any(is.infinite(unlist(p)))) warning("Proportions have NA or Inf - check results!",call. = FALSE)
      if(any(is.na(r))|any(is.infinite(unlist(r)))) warning("Transformed correlations have NA or Inf - check results!",call. = FALSE)
      if(any(is.na(logOR))|any(is.infinite(unlist(logOR)))) warning("Log odds ratios have NA or Inf - check results!",call. = FALSE)
    }
    fitcheck(est.result())
  })

  output$itfit <- shiny::renderPrint({
    if (input$goButton == 0)
      return()
    itf()
  })

  itfplot <- shiny::reactive({
    itemfit(est.result())
  })


  ip <- shiny::reactive({
    if (input$goButton == 0) return()
    coef(est.result(),what = input$ips,withSE=input$ipse)
  })

  output$ip <- shiny::renderPrint({
    if (input$goButton == 0)
      return()
    coef(est.result(),what = input$ips,withSE=input$ipse)
  })

  output$pparm <- shiny::renderPrint({
    head(personparm(object = est.result(),what = input$pp),10)
  })

  output$plc.output <- shiny::renderPrint({
    x=extract(est.result(),"posterior.prob")
    xx <- data.frame(latentclass=attr(x,"dimnames")[[2]],proportion=c(x))
    if(input$plc=="default"){
      return(head(xx,10))
    }else if(input$plc=="decreasing"){
      return(head(xx[order(xx$proportion,decreasing = TRUE),],10))
    }else if(input$plc=="increasing"){
      return(head(xx[order(xx$proportion,decreasing = FALSE),],10))
    }

  })

  q <- shiny::reactive({
    if (input$qvalcheck == 0)  return()
    GDINA::Qval(est.result(),method = input$qv.method,eps = input$PVAFcutoff)
  })
  output$sugQ <- shiny::renderPrint({
    if (input$qvalcheck == 0)  return()
    # extract(q(),what = "sug.Q")
    q()
  })

  m <- shiny::reactive({
    if (input$modelsel == 0)  return()
    modelcomp(est.result())
  })

  output$ws <- shiny::renderPrint({
    if (input$modelsel == 0)  return()
    extract(m(),what = "stats")
  })
  output$pv <- shiny::renderPrint({
    if (input$modelsel == 0)  return()
    extract(m(),what = "pvalues")
  })
  output$ws <- shiny::renderPrint({
    if (input$modelsel == 0)  return()
    extract(m(),what = "stats")
  })
  output$ss <- shiny::renderPrint({
    if (input$modelsel == 0)  return()
    print(m())
  })


  #########################
  # IRF plots
  #
  #########################


makeIRFplot <- function(){
  if (input$goButton == 0)
    return()

  if (input$item.plot<=extract(est.result(),"ncat"))
    plot(est.result(),item = input$item.plot, withSE = input$IRFplotse)
}

output$plot <- shiny::renderPlot({
  if (input$goButton == 0)
    return()
  makeIRFplot()
})

output$downloadIRFplot <- shiny::downloadHandler(
  filename = function() {
    paste('IRFPlot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    print(makeIRFplot())
    dev.off()
  }
)

#########################
# individual mastery plots
#
#########################


makeMpplot <- function(){
  if (input$goButton == 0)
    return()
  df <- personparm(est.result(),"mp")
  att.names <- colnames(df)
  person <- as.numeric(unlist(strsplit(input$personid,",")))
  np <- length(person)
  if(np>1){
    dff <- c(t(df[person,]))
    dat <- data.frame(att = rep(att.names,np),mp = dff,person = factor(rep(person,each = ncol(df))))
    x <- ggplot2::ggplot(data = dat, ggplot2::aes_string(x = "att", y = "mp")) +
      ggplot2::geom_bar(stat = "identity", position = "dodge",ggplot2::aes_string(fill = "person")) +
      ggplot2::ylim(0,1)+
      ggplot2::labs(x = "Attributes", y = "Mastery probabilities",
                 title = paste("Mastery probabilities"))
  }else{
    dff <- c(df[person,])
    dat <- data.frame(att = att.names,mp = dff,person = factor(rep(person,ncol(df))))
    x <- ggplot2::ggplot(data = dat, ggplot2::aes_string(x = "att", y = "mp")) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::ylim(0,1)+
      ggplot2::labs(x = "Attributes", y = "Mastery probabilities",
                    title = paste("Mastery probabilities for individual",person))
  }
  if(input$HPlot){
    print(x + ggplot2::coord_flip())
  }else{
    print(x)
  }
}

output$Mplot <- shiny::renderPlot({
  if (input$goButton == 0)
    return()
  makeMpplot()
})

output$downloadmpplot <- shiny::downloadHandler(
  filename = function() {
    paste('MPPlot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    makeMpplot()
    dev.off()
  }
)

#########################
# individual posterior plots
#
#########################

makeiPostProbplot <- function(){
  x <- exp(indlogPost(est.result()))[input$ippid,]
  lc.names <- attr(x,"names")
  if (input$ippid>extract(est.result(),"nobs"))
    return(NULL)
  nc <- min(input$inlc,length(x))
  xx <- data.frame(LC = lc.names,prob = c(x))
  if(input$ippplc=="default"){
    y <- xx[seq_len(nc),]
    # y <- y[complete.cases(y),]
    if(input$ippAdir){
      print(ggplot2::ggplot(data=y, ggplot2::aes(x=LC, y=prob)) +
              ggplot2::geom_bar(stat="identity")+
              ggplot2::coord_flip()+ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
    }else{
      print(ggplot2::ggplot(data=y, ggplot2::aes(x=LC, y=prob)) +
              ggplot2::geom_bar(stat="identity")+
              ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
    }
  }else if(input$ippplc=="decreasing"){
    y <- xx[order(-c(x))[seq_len(nc)],]

    if(input$ippAdir){
      print(ggplot2::ggplot(data=y, ggplot2::aes(x=reorder(LC,-prob), y=prob)) +
              ggplot2::geom_bar(stat="identity")+
              ggplot2::coord_flip()+ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
    }else{
      print(ggplot2::ggplot(data=y, ggplot2::aes(x=reorder(LC,-prob), y=prob)) +
              ggplot2::geom_bar(stat="identity")+
              ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
    }
  }else if(input$ippplc=="increasing"){
    y <- xx[order(c(x))[seq_len(nc)],]

    if(input$ippAdir){
      print(ggplot2::ggplot(data=y, ggplot2::aes(x=reorder(LC,prob), y=prob)) +
              ggplot2::geom_bar(stat="identity")+
              ggplot2::coord_flip()+ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
    }else{
      print(ggplot2::ggplot(data=y, ggplot2::aes(x=reorder(LC,prob), y=prob)) +
              ggplot2::geom_bar(stat="identity")+
              ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
    }
  }


}
output$iPostProbplot <- shiny::renderPlot({
  if (input$goButton == 0)
    return()
  makeiPostProbplot()
})

output$downloadiPPplot <- shiny::downloadHandler(
  filename = function() {
    paste('iPP-Plot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    makeiPostProbplot()
    dev.off()
  }
)
#########################
# Group level
# attribute prevalence plots
#
#########################

makeAPpplot <- function(){
  if (input$goButton == 0)
    return()
  x <- extract(est.result(),"prevalence")[[1]]

  l <- rev(paste0("Level",seq_len(ncol(x))-1))

  df <- data.frame(Attribute=rep(rownames(x),ncol(x)),
                   Levels=factor(rep(colnames(x),each=nrow(x)),labels = l,levels = l,ordered = TRUE),
                   Prevalence=round(c(x),3),
                   label_ypos=c(t(apply(x,1,cumsum))))

  if(input$Adir){
    print(ggplot2::ggplot(data=df, ggplot2::aes(x=Attribute, y=Prevalence, fill=Levels)) +
            ggplot2::geom_bar(stat="identity")+
            ggplot2::geom_text(ggplot2::aes(y=label_ypos, label=Prevalence), hjust = 1.6,
                      color="white", size=4.5)+
            ggplot2::scale_fill_brewer(palette=input$palette)+
            ggplot2::theme_minimal() + ggplot2::coord_flip())
  }else{
    print(ggplot2::ggplot(data=df, ggplot2::aes(x=Attribute, y=Prevalence, fill=Levels)) +
            ggplot2::geom_bar(stat="identity")+
            ggplot2::geom_text(ggplot2::aes(y=label_ypos, label=Prevalence), vjust=1.6,
                      color="white", size=4.5)+
            ggplot2::scale_fill_brewer(palette=input$palette)+
            ggplot2::theme_minimal())
  }
}

output$APplot <- shiny::renderPlot({
  if (input$goButton == 0)
    return()
  makeAPpplot()
})

output$downloadAPplot <- shiny::downloadHandler(
  filename = function() {
    paste('AP-Plot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    makeAPpplot()
    dev.off()
  }
)



######################
#
# Group probability of each latent class
#
######################

makePostProbplot <- function(){
  x <- extract(est.result(),"posterior.prob")
  xx <- data.frame(LC=c(attr(x,"dimnames")[[2]]),prob=round(c(x),4))
  nc <- min(input$nlc,nrow(xx))
    if(input$ppplc=="default"){
      y <- xx[seq_len(nc),]
      y <- y[complete.cases(y),]
      if(input$ppAdir){
        print(ggplot2::ggplot(data=y, ggplot2::aes(x=LC, y=prob)) +
                ggplot2::geom_bar(stat="identity")+
                ggplot2::coord_flip()+ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
      }else{
        print(ggplot2::ggplot(data=y, ggplot2::aes(x=LC, y=prob)) +
                ggplot2::geom_bar(stat="identity")+
                ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
      }
    }else if(input$ppplc=="decreasing"){
      y <- xx[order(-x)[seq_len(nc)],]
      y <- y[complete.cases(y),]
      if(input$ppAdir){
        print(ggplot2::ggplot(data=y, ggplot2::aes(x=reorder(LC,-prob), y=prob)) +
                ggplot2::geom_bar(stat="identity")+
                ggplot2::coord_flip()+ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
      }else{
        print(ggplot2::ggplot(data=y, ggplot2::aes(x=reorder(LC,-prob), y=prob)) +
                ggplot2::geom_bar(stat="identity")+
                ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
      }
    }else if(input$ppplc=="increasing"){
      y <- xx[order(x)[seq_len(nc)],]
      y <- y[complete.cases(y),]
      if(input$ppAdir){
        print(ggplot2::ggplot(data=y, ggplot2::aes(x=reorder(LC,prob), y=prob)) +
                ggplot2::geom_bar(stat="identity")+
                ggplot2::coord_flip()+ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
      }else{
        print(ggplot2::ggplot(data=y, ggplot2::aes(x=reorder(LC,prob), y=prob)) +
                ggplot2::geom_bar(stat="identity")+
                ggplot2::ylab("Posterior probability")+ggplot2::xlab("Latent Class"))
      }
    }


}

output$PostProbplot <- shiny::renderPlot({
  if (input$goButton == 0)
    return()
  makePostProbplot()
})

output$downloadPPplot <- shiny::downloadHandler(
  filename = function() {
    paste('PP-Plot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    makePostProbplot()
    dev.off()
  }
)


#########################
# mesa plots
#
#########################


makeMesaplot <- function(){
  if (input$qvalcheck == 0)  return()
  plot(q(),item = input$item.mesaplot,type = input$mesatype, data.label = input$datalabel)
}

output$mesaplot <- shiny::renderPlot({
    if (input$qvalcheck == 0)  return()
    makeMesaplot()
  })
output$downloadMesaplot <- shiny::downloadHandler(
  filename = function() {
    paste('MesaPlot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    print(makeMesaplot())
    dev.off()
  }
)

#########################
# heat plots
#
#########################

makeHeatplot <- function(){

    item.pair.1 <- item.pair.2 <- unadj.pvalue <- test.adj.pvalue <- NULL
    if(input$heatmap.type=="log odds ratio"){
      df <- extract(itfplot(),"logOR")
    }else{
      df <- extract(itfplot(),"r")
    }

    if(input$heatmap.adjust){
      p <- ggplot2::ggplot(df, ggplot2::aes(x=factor(item.pair.2),
                                   y=factor(item.pair.1),
                                   fill=test.adj.pvalue))+
        ggplot2::geom_tile()+ ggplot2::scale_fill_gradient(low="red",
                                         high="gray",
                                         limits=c(0,0.05))+
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Items", y = "Items",
                      title = paste("Heatmap plot for adjusted p-values of ",input$heatmap.type))
    }else{
      p <- ggplot2::ggplot(df, ggplot2::aes(x=factor(item.pair.2),
                                   y=factor(item.pair.1),
                                   fill=unadj.pvalue))+
        ggplot2::geom_tile()+ ggplot2::scale_fill_gradient(low="red",
                                         high="gray",
                                         limits=c(0,0.05))+
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Items", y = "Items",
             title = paste("Heatmap plot for unadjusted p-values of ",input$heatmap.type))
    }


    print(p)
  }
output$heatplot <- shiny::renderPlot({
    if (input$goButton == 0)
      return()
    makeHeatplot()
  })

output$downloadHeatPlot <- shiny::downloadHandler(
  filename = function() {
    paste('HeatPlot', Sys.Date(), '.pdf', sep='')
  },
  content = function(FILE=NULL) {
    pdf(file=FILE)
    makeHeatplot()
    dev.off()
  }
)


  output$downloadpp <- shiny::downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      ext <- switch(input$ppfiletype, "csv" = ".csv", "tsv" = ".txt")
      paste(input$pp, ext, sep = "")

    },

    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      sep <- switch(input$ppfiletype, "csv" = ",", "tsv" = "\t")

      write.table(personparm(object = est.result(),what = input$pp), file, sep = sep,
                   row.names = FALSE)
    }
  )
  output$downloadplc <- shiny::downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      ext <- switch(input$plcfiletype, "csv" = ".csv", "tsv" = ".txt")
      paste(input$plc, ext, sep = "")
    },

    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      sep <- switch(input$plcfiletype, "csv" = ",", "tsv" = "\t")
      x=extract(est.result(),"posterior.prob")
      xx <- data.frame(latentclass=attr(x,"dimnames")[[2]],proportion=c(x))
      if(input$plc=="default"){
        y <- xx
      }else if(input$plc=="decreasing"){
        y <- xx[order(xx$proportion,decreasing = TRUE),]
      }else if(input$plc=="increasing"){
        y <- xx[order(xx$proportion,decreasing = FALSE),]
      }
      # Write to a file specified by the 'file' argument
      write.table(y, file, sep = sep, row.names = FALSE)
    }
  )




  })


