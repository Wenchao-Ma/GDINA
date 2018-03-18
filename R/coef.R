#' @include GDINA.R
#' @export
#' @description To calculate structural parameters for item and joint attribute distributions, use method \code{\link{coef}}.
#' @describeIn GDINA extract structural parameter estimates
coef.GDINA <-
  function(object,
           what = c("catprob","delta","gs","itemprob","LCprob","rrum","lambda"),
           withSE = FALSE, SE.type = 2,digits = 4, ...)
  {

    if(!class(object)=="GDINA") stop("object must be a GDINA estimate.",call. = FALSE)
    what <- match.arg(what)
    if(extract(object,"att.str")){
      if(tolower(what)=="delta")stop("Delta parameters are not availabel for models with structured attributes.",call. = FALSE)
    }
    if(tolower(what)=="catprob"){
      if(withSE){
        out <- mapply(rbind,extract(object,what = "catprob.parm"),
                      extract(object,what = "catprob.se", SE.type = SE.type),SIMPLIFY = F)
        out <- lapply(out,function(x){rownames(x) <- c("Est.","S.E.");round(x,digits)})
      }else{
        out <- lapply(extract(object,what = "catprob.parm"),round,digits)
      }
    }else if(tolower(what)=="itemprob"){
      if(withSE&!extract(object,what = "sequential")){
        out <- mapply(rbind,extract(object,what = "itemprob.parm"),
                      extract(object,what = "itemprob.se", SE.type = SE.type),SIMPLIFY = F)
        out <- lapply(out,function(x){rownames(x) <- c("Est.","S.E.");round(x,digits)})
      }else{
        out <- lapply(extract(object,what = "itemprob.parm"),round,digits)
      }
    }else if(tolower(what)=="lcprob"){

      out <- round(extract(object,what = "LCprob.parm"),digits)

    }else if(tolower(what)=="gs"){
      if(withSE){
        gs <- t(sapply(extract(object,what = "catprob.parm"),
                       function(x)c(x[1],1-rev(x)[1])))
        se <- t(sapply(extract(object,what = "catprob.se", SE.type = SE.type),
                       function(x)c(x[c(1,length(x))])))
        out <- round(cbind(gs,se),digits)
        colnames(out) <- c("guessing","slip","SE[guessing]","SE[slip]")
      }else{
        out <- round(t(sapply(extract(object,what = "catprob.parm"),
                              function(x)c(x[1],1-rev(x)[1]))),digits)
        colnames(out) <- c("guessing","slip")
      }
    }else if(tolower(what)=="delta"){
      d <- format_delta(delta = extract(object,what = "delta.parm"),
                        model = extract(object,what = "models_numeric"),
                        Kj = extract(object,what = "Kj"),
                        item.names = extract(object,what = "item.names"),
                        digits = digits+1)
      if(withSE){
        out <- mapply(rbind,d,
                      extract(object,what = "delta.se", SE.type = SE.type),SIMPLIFY = F)
        out <- lapply(out,function(x){rownames(x) <- c("Est.","S.E.");round(x,digits)})
      }else{
        out <- lapply(d,round,digits)
      }
    }else if(tolower(what)=="rrum"){
      J <- extract(object,what = "ncat")
      ip <- extract(object,what = "catprob.parm")
      models <- extract(object,what = "models_numeric")
      if(all(models!=5)) stop("RRUM parameters are only available when RRUM is fitted.",call. = FALSE)
      out <- vector("list",J)
      for (j in 1:J){
        if (models[j]==5){
          parj <- ip[[j]]
          Lj <- length(parj)
          Kj <- log(Lj,2)
          pi.star <- parj[Lj]
          out[[j]] <- round(c(pi.star,rev(parj[(Lj-Kj):(Lj-1)])/pi.star),digits)
          names(out[[j]]) <- c("pi*",paste0("r",1:Kj))
        }

      }
      names(out) <- extract(object,what = "item.names")
      if(withSE) message("Standard errors are not available for RRUM parameters.")
    }else if(tolower(what)=="lambda"){
      if(extract(object,"ngroup")==1){
        out <- round(extract(object,"struc.parm")[[1]],digits)
        if(any(extract(object,"att.dist")=="higher.order")){
          rownames(out) <- paste0("A",seq_len(extract(object,"natt")))
          colnames(out) <- c("slope","intercept")
        }else if(extract(object,"att.dist")=="saturated"){
          lab <- paste0("p(",apply(attributepattern(Q=extract(object,"Q")),1,paste0,collapse=""),")")
          names(out) <- lab
        }else if(extract(object,"att.dist")=="independent"){
          names(out) <- paste0("P(A",seq_len(extract(object,"natt")),")")
        }
      }else{
        out <- lapply(extract(object,"struc.parm"),round,digits=digits)
        names(out) <- paste0("G",seq_len(extract(object,"ngroup")))
        for(g in seq_len(extract(object,"ngroup"))){
          if(extract(object,"att.dist")[g]=="higher.order"){
            rownames(out[[g]]) <- paste0("A",seq_len(extract(object,"natt")))
            colnames(out[[g]]) <- c("slope","intercept")
          }else if(extract(object,"att.dist")[g]=="saturated"){
            lab <- paste0("p(",apply(attributepattern(Q=extract(object,"Q")),1,paste0,collapse=""),")")
            names(out[[g]]) <- lab
          }else if(extract(object,"att.dist")[g]=="independent"){
            names(out[[g]]) <- paste0("P(A",seq_len(extract(object,"natt")),")")
          }
        }


      }
      if(withSE) warning("Please Calculate SEs for parameters of the structural model using bootstrap method.",call. = FALSE)
    }else{
      stop(sprintf("No structural parameters called \'%s\'", what), call.=FALSE)
    }

    return(out)
  }
