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

  est.result <- shiny::reactive(
    withProgress(message = 'Model Estimating', value = 0.5, {
    inFile1 <- input$file1
    dat <- read.csv(inFile1$datapath, header = input$header,
             sep = input$sep, quote = input$quote)

    inFile2 <- input$file2
    Q <- read.csv(inFile2$datapath, header = input$header2,
                    sep = input$sep2, quote = input$quote)
    if(input$attdis==0){
      HOdist <- "saturated"
    }else if(input$attdis==1){
      HOdist <- "higher.order"
    }else if(input$attdis==2){
      HOdist <- "fixed"
    }

      est <- GDINA::GDINA(dat = dat, Q = Q, model = input$type,
                          verbose = 0,att.dist = HOdist,
                          higher.order = list(model = input$hom),
                          sequential = input$seq,
                          mono.constraint = input$mono)


    est
  }))


  ##################
  # Summary
  ##################
  info <- shiny::reactive({
    summary(est.result())
  })
  iter.info <- shiny::reactive({
    est.info <- function(x) {
      cat("\nThe Generalized DINA Model Framework  \n")
      packageinfo <- utils::packageDescription("GDINA")
      cat( paste( "   GDINA Version " , packageinfo$Version , " (" , packageinfo$Date , ")" , sep="") , "\n" )
      cat(  "   Wenchao Ma & Jimmy de la Torre \n" )

      cat("\nNumber of items       =", extract(x,"nitem"), "\n")
      cat("Number of individuals =", extract(x,"nobs"), "\n")
      cat("Number of attributes  =", extract(x,"natt"), "\n")
      M <- c("GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM")
      cat("Number of iterations  =", extract(x,"nitr"), "\n")
      cat("Fitted model(s)       =\n")
      print(extract(x,"models"))

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




  output$info <- shiny::renderPrint({
    if (input$goButton == 0)
      return()
    info()
  })

  output$iter.info <- shiny::renderPrint({
    if (input$goButton == 0)
      return()
    iter.info()
  })


  itf <- shiny::reactive({
    itemfit(est.result())
  })

  output$itfit <- shiny::renderPrint({
    print(itf())
  })







  ip <- shiny::reactive({
    if (input$goButton == 0) return()
    coef(est.result(),what = input$ips,withSE=TRUE)
  })

  output$ip <- shiny::renderPrint({
    if (input$goButton == 0)
      return()
    coef(est.result(),what = input$ips,withSE=TRUE)
  })

  output$pparm <- shiny::renderPrint({
    head(personparm(object = est.result(),what = input$pp),10)
  })

  q <- shiny::reactive({
    if (input$qvalcheck == 0)  return()
    Qval(est.result(),eps = input$PVAFcutoff)
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

makeIRFplot <- function(){
  if (input$goButton == 0)
    return()
  inFile2 <- input$file2
  Q <- read.csv(inFile2$datapath, header = input$header,
                sep = input$sep, quote = input$quote)
  if (input$item.plot<1||input$item.plot>nrow(Q)) NULL
  plot(est.result(),IRF.args = list(item = input$item.plot, errorbar = input$IRFplotse))
}

output$plot <- shiny::renderPlot({
  if (input$goButton == 0)
    return()
  makeIRFplot()
})

makeMesaplot <- function(){
  if (input$qvalcheck == 0)  return()
  plot(q(),item = input$item.mesaplot,type = input$mesatype, data.label = input$datalabel)
}

  output$mesaplot <- shiny::renderPlot({
    if (input$qvalcheck == 0)  return()
    makeMesaplot()
  })

  makeHeatplot <- function(){

    item.pair.1 <- item.pair.2 <- unadj.pvalue <- test.adj.pvalue <- NULL
    if(input$heatmap.type=="log odds ratio"){
      df <- extract(itf(),"logOR")
    }else{
      df <- extract(itf(),"r")
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
  output$heatplot1 <- shiny::renderPlot({
    if (input$goButton == 0)
      return()
    makeHeatplot()
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


  output$downloadpp <- shiny::downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste(input$pp, input$ppfiletype, sep = ".")
    },

    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      sep <- switch(input$ppfiletype, "csv" = ",", "tsv" = "\t")

      # Write to a file specified by the 'file' argument
      write.table(personparm(object = est.result(),what = input$pp), file, sep = sep,
                  row.names = FALSE)
    }
  )
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



  })


