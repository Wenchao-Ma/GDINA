#' Calculating standard errors and variance-covariance matrix using bootstrap methods
#'
#' This function conducts nonparametric and parametric bootstrap to calculate standard errors of model parameters.
#' Parametric bootstrap is only applicable for single group models.
#'
#' @param GDINA.obj an object of class GDINA
#' @param bootsample the number of bootstrap samples
#' @param type type of bootstrap method. Can be \code{parametric} or \code{nonparametric}
#' @param randomseed random seed for resampling
#'
#' @return itemparm.se standard errors for item probability of success in list format
#' @return delta.se standard errors for delta parameters in list format
#' @return strucparm.se standard errors for structural model parameters
#' @return boot.est resample estimates
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @export
#' @examples
#' \dontrun{
#' # For illustration, only 5 resamples are run
#' # results are definitely not reliable
#'
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA",att.dist = "higher.order")
#' boot.fit <- bootSE(fit,bootsample = 5,randomseed=123)
#' boot.fit$delta.se
#' }
bootSE <- function(GDINA.obj,bootsample=50,type = "nonparametric",randomseed=12345){

  if (exists(".Random.seed", .GlobalEnv)) oldseed <- .GlobalEnv$.Random.seed else  oldseed <- NULL

  set.seed(randomseed)

  Y <- extract(GDINA.obj,"dat")
  Q <- extract(GDINA.obj,"Q")
  N <- nrow(Y)
  J <- ncol(Y)
  K <- ncol(Q)
  no.mg <- extract(GDINA.obj,"ngroup")
  out.list <- vector("list",bootsample)
  # By default

  GDINA.options <- formals(GDINA)
  GDINA.options <- GDINA.options[-c(1,2,length(GDINA.options))]
  # user specified
  tmp <- as.list(GDINA.obj$extra$call)[-c(1:3)]
  GDINA.options[names(GDINA.options)%in%names(tmp)] <- tmp
  GDINA.options$verbose <- 0

    for (r in 1:bootsample){
      if(tolower(type)=="parametric"){
        simdat <- simGDINA(N, Q, catprob.parm = GDINA.obj$catprob.parm,
                           attribute = attributepattern(Q = Q)[sample(1:2^K,N,replace = TRUE,prob = GDINA.obj$posterior.prob),])$dat

      }else if(tolower(type)=="nonparametric"){
        simdat <- Y[sample(1:N,N,replace = TRUE),]
      }

      boot.out <- do.call(GDINA,c(list(dat = simdat,Q = Q),GDINA.options))
      out.list[[r]] <- list(itemprob = boot.out$catprob.parm,delta = boot.out$delta.parm,
                            jointAtt = boot.out$struc.parm)
    }


  ip <- out.list[[1]]$itemprob
  d <- out.list[[1]]$delta
  for(r in 2:(bootsample)){
    ip <- Map(rbind,ip,out.list[[r]]$itemprob)
    d <- Map(rbind,d,out.list[[r]]$delta)
  }
  se.ip <- lapply(ip,function(x)apply(x,2,sd))
  se.d <- lapply(d,function(x)apply(x,2,sd))
  se.jointAtt <- list()
  for(g in 1:no.mg){
    tmp <- c(out.list[[1]]$jointAtt[[g]])
    for(r in 2:(bootsample)){
      tmp <- rbind(tmp,c(out.list[[r]]$jointAtt[[g]]))
    }
    se.tmp <- apply(tmp,2,sd)
    if(extract(GDINA.obj,"att.dist")[g]=="higher.order"){
      se.jointAtt[[g]] <- matrix(se.tmp,ncol = 2)
    }else{
      se.jointAtt[[g]] <- se.tmp
    }
  }



  if (!is.null(oldseed)) .GlobalEnv$.Random.seed <- oldseed  else  rm(".Random.seed", envir = .GlobalEnv)

  return(list(itemparm.se=se.ip,delta.se=se.d,strucparm.se=se.jointAtt,boot.est=out.list))
}

