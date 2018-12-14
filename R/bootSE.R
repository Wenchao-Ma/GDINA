#' Calculating standard errors and variance-covariance matrix using bootstrap methods
#'
#' This function conducts nonparametric and parametric bootstrap to calculate standard errors of model parameters.
#' Parametric bootstrap is only applicable to single group models.
#'
#' @param GDINA.obj an object of class GDINA
#' @param bootsample the number of bootstrap samples
#' @param type type of bootstrap method. Can be \code{parametric} or \code{nonparametric}
#' @param randomseed random seed for resampling
#'
#' @return itemparm.se standard errors for item probability of success in list format
#' @return delta.se standard errors for delta parameters in list format
#' @return lambda.se standard errors for structural parameters of joint attribute distribution
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
#' boot.fit$lambda.se
#' }
bootSE <- function(GDINA.obj,bootsample=50,type = "nonparametric",randomseed=12345){

  if (exists(".Random.seed", .GlobalEnv)) oldseed <- .GlobalEnv$.Random.seed else  oldseed <- NULL

  if(!is.null(extract(GDINA.obj,"att.str")))
    stop("bootSE is not available for models with attribute structures.",call. = FALSE)
  set.seed(randomseed)

  Y <- extract(GDINA.obj,"dat")
  Q <- extract(GDINA.obj,"Q")
  N <- nrow(Y)
  J <- ncol(Y)
  K <- ncol(Q)
  no.mg <- extract(GDINA.obj,"ngroup")
  stopifnot(no.mg==1)
  lambda <- delta <- itemprob <- vector("list",bootsample)
  # By default

  GDINA.options <- formals(GDINA)
  GDINA.options <- GDINA.options[-c(1,2,length(GDINA.options))]
  # user specified
  tmp <- as.list(GDINA.obj$extra$call)[-c(1:3)]
  GDINA.options[names(GDINA.options)%in%names(tmp)] <- tmp
  GDINA.options$verbose <- 0
  att <- extract(GDINA.obj,"attributepattern")
    for (r in 1:bootsample){

      if(tolower(type)=="parametric"){

        simdat <- simGDINA(N, Q, catprob.parm = GDINA.obj$catprob.parm,
                           attribute = att[sample(seq_len(nrow(att)),N,replace = TRUE,prob = GDINA.obj$posterior.prob),])$dat

      }else if(tolower(type)=="nonparametric"){
        simdat <- Y[sample(1:N,N,replace = TRUE),]
      }

      boot.out <- do.call(GDINA,c(list(dat = simdat,Q = Q),GDINA.options))
      lambda[[r]] <- c(boot.out$struc.parm)
      itemprob[[r]]<- boot.out$catprob.parm
      delta[[r]]<- boot.out$delta.parm
    }



  se.ip <- lapply(do.call(Map,c(f="rbind",itemprob)),function(x)apply(x,2,sd))
  se.d <- lapply(do.call(Map,c(f="rbind",delta)),function(x)apply(x,2,sd))
  se.lambda <- apply(do.call(rbind,lambda),2,sd)

    if(extract(GDINA.obj,"att.dist")=="higher.order"){
      se.jointAtt <- matrix(se.lambda,ncol = 2)
    }else{
      se.jointAtt <- se.lambda
    }




  if (!is.null(oldseed)) .GlobalEnv$.Random.seed <- oldseed  else  rm(".Random.seed", envir = .GlobalEnv)

  return(list(itemparm.se=se.ip,delta.se=se.d,lambda.se=se.jointAtt,boot.est=list(lambda=lambda,itemprob=itemprob,delta=delta)))
}

