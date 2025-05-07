#' Calculating standard errors and variance-covariance matrix using bootstrap methods
#'
#' This function conducts nonparametric and parametric bootstrap to calculate standard errors of model parameters.
#' Parametric bootstrap is only applicable to single group models.
#'
#' @param GDINA.obj an object of class GDINA
#' @param bootsample the number of bootstrap samples
#' @param type type of bootstrap method. Can be \code{parametric} or \code{nonparametric}
#' @param cores number of cores to be used for calculation; Default = 1; if \code{cores>1}, bootstrap will run in parallel with \code{foreach} package.
#' @param seed random seed for resampling
#'
#' @return itemparm.se standard errors for item probability of success in list format
#' @return delta.se standard errors for delta parameters in list format
#' @return lambda.se standard errors for structural parameters of joint attribute distribution
#' @return boot.est resample estimates
#'
#' @importFrom foreach %dopar% foreach
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#' @references
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive Diagnosis Modeling. \emph{Journal of Statistical Software, 93(14)}, 1-26.
#'
#' @export
#' @examples
#' \dontrun{
#' # For illustration, only 5 resamples are run
#' # results are definitely not reliable
#'
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA",att.dist = "higher.order")
#' boot.fit <- bootSE(fit,bootsample = 5,seed=123)
#' boot.fit$delta.se
#' boot.fit$lambda.se
#' }
bootSE <- function(GDINA.obj,bootsample=50,type = "nonparametric",cores = 1, seed=12345){

  if (exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- oldseed)
  }else{
    on.exit(rm(".Random.seed", envir = .GlobalEnv))
  }

  stopifnot(isa(GDINA.obj,"GDINA"))
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
  if(cores == 1){
    set.seed(seed)
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

  }else{
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("doParallel is required for parallel execution.")
    }
    if (!requireNamespace("doRNG", quietly = TRUE)) {
      stop("doRNG is required for reproducibility in parallel execution.")
    }

    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(seed)

    fe <- foreach(r = 1 : bootsample,.packages = "GDINA") %dopar% {

              if(tolower(type)=="parametric"){

                simdat <- simGDINA(N, Q, catprob.parm = GDINA.obj$catprob.parm,
                                   attribute = att[sample(seq_len(nrow(att)),N,replace = TRUE,prob = GDINA.obj$posterior.prob),])$dat

              }else if(tolower(type)=="nonparametric"){
                simdat <- Y[sample(1:N,N,replace = TRUE),]
              }

              boot.out <- do.call(GDINA,c(list(dat = simdat,Q = Q),GDINA.options))
              list(lambda = c(boot.out$struc.parm),
              itemprob = boot.out$catprob.parm,
              delta = boot.out$delta.parm)

    }
    lambda <- lapply(fe,function(x)x[[1]])
    itemprob <- lapply(fe,function(x)x[[2]])
    delta <- lapply(fe,function(x)x[[3]])
    parallel::stopCluster(cl)
  }



  se.ip <- lapply(do.call(Map,c(f="rbind",itemprob)),function(x)apply(x,2,sd))
  se.d <- lapply(do.call(Map,c(f="rbind",delta)),function(x)apply(x,2,sd))
  se.lambda <- apply(do.call(rbind,lambda),2,sd)

    if(extract(GDINA.obj,"att.dist")=="higher.order"){
      se.jointAtt <- matrix(se.lambda,ncol = 2)
    }else{
      se.jointAtt <- se.lambda
    }


return(list(itemparm.se=se.ip,delta.se=se.d,lambda.se=se.jointAtt,boot.est=list(lambda=lambda,itemprob=itemprob,delta=delta)))
}

