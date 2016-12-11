#' Calculating standard errors and variance-covariance matrix using parametric bootstrap method
#'
#' This function conducts parametric bootstrap to calculate standard errors and variance-covariance matrix
#' for item probability of success, delta parameters and higher order parameters (if necessary)
#'
#' @param GDINA.obj an object of class GDINA
#' @param Rsample the number of bootstrap samples
#' @param digits How many decimal places in each number? The default is 4.
#'
#' @return item.prob.se standard errors for item probability of success in list format
#' @return delta.se standard errors for delta parameters in list format
#' @return HO.slope.se standard errors for slope parameters of higher order model
#' @return HO.intercept.se standard errors for intercept parameters of higher order model
#' @return item.prob.cov variance-covariance matrix for item probability of success
#' @return delta.cov variance-covariance matrix for delta parameters
#' @return HO.slope.var variance-covariance matrix for slope
#' @return HO.intercept.var  variance-covariance matrix for intercept
#' @return item.prob.resamples resample estimates for item probability of success
#' @return delta.resamples resample estimates for delta parameters
#' @return HO.slope.resamples resample estimates for slope
#' @return HO.intercept.resamples resample estimates for intercept
#'
#' @export
#' @examples
#'
#' # For illustration, only 5 resamples are run
#' # results are definitely not reliable
#'
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA",higher.order = TRUE)
#' boot.mod1 <- bootSE(mod1,Rsample = 50,parallel=TRUE)
#' boot.mod1$item.prob.se
#' boot.mod1$HO.intercept.se
#'
#' mod2 <- GDINA(dat = dat, Q = Q, model = "GDINA",higher.order = FALSE)
#' boot.mod2 <- bootSE(mod2,Rsample = 200)
#' boot.mod2$item.prob.se
#'
bootSE <- function(GDINA.obj,Rsample=100,digits = 4,
                   conv.crit = 0.001, type = "parametric",
                   parallel=FALSE){
  Y <- GDINA.obj$options$dat
  Q <- GDINA.obj$options$Q
  N <- nrow(Y)
  J <- ncol(Y)
  K <- ncol(Q)

  out.list <- vector("list",Rsample)
  GDINA.options <- list(model = GDINA.obj$options$model,
                        higher.order = GDINA.obj$options$higher.order,
                        higher.order.model = GDINA.obj$options$higher.order.model,
                        higher.order.method = GDINA.obj$options$higher.order.method,
                        catprob.matrix = GDINA.obj$catprob.matrix,
                        mono.constraint = GDINA.obj$options$mono.constraint,
                        empirical = GDINA.obj$options$empirical,
                        att.prior = GDINA.obj$options$att.prior,
                        att.str = GDINA.obj$options$att.str,
                        conv.crit = conv.crit,
                        maxitr = 1000,
                        higher.order.struc.parm = GDINA.obj$options$higher.order.att.parm)
if(parallel){
  if(require(foreach)){
    cl <- parallel::makeCluster(parallel::detectCores())
    doParallel::registerDoParallel(cl)

    out.list <- foreach::foreach (r=1:Rsample,.packages = "GDINA")%dopar%{
      if(tolower(type)=="parametric"){
         simdat <- simGDINA(N, Q, catprob.parm = GDINA.obj$catprob.parm,
                          attribute = attributepattern(K,T,Q)[sample(1:2^K,N,replace = TRUE,prob = GDINA.obj$posterior.prob),])$dat
        # simdat <- GDINA.sim(N, Q, itemprob.param = GDINA.obj$itemprob.param)$dat
      }else if(tolower(type)=="nonparametric"){
        simdat <- Y[sample(1:N,N,replace = TRUE),]
      }
      boot.out <- GDINA(simdat,Q,model = GDINA.options$model, higher.order = GDINA.options$higher.order,
                        higher.order.model = GDINA.options$higher.order.model,higher.order.method = GDINA.options$higher.order.method,
                        verbose = FALSE, catprob.parm = GDINA.obj$catprob.parm,
                        mono.constraint = GDINA.options$mono.constraint,
                        empirical = GDINA.options$empirical, att.prior = GDINA.options$att.prior,
                        att.str = GDINA.options$att.str,
                        nstarts = 3, conv.crit = GDINA.options$conv.crit,
                        maxitr = GDINA.options$maxitr, higher.order.struc.parm =  GDINA.options$higher.order.parm)
      return(list(itemprob = unlist(boot.out$catprob.parm),delta = unlist(boot.out$delta.parm), HO = c(boot.out$higher.order.struc.parm)))
    }
  }else{
    stop("foreach package is required if parallel=TRUE.",call. = FALSE)
  }
parallel::stopCluster(cl)
}else{
  for (r in 1:Rsample){
    sim <- simGDINA(N, Q, catprob.parm = GDINA.obj$catprob.parm,
                     attribute = attributepattern(K,T,Q)[sample(1:2^K,N,replace = TRUE,prob = GDINA.obj$posterior.prob),])
    boot.out <- GDINA(sim$dat,Q,model = GDINA.options$model, higher.order = GDINA.options$higher.order,
                      higher.order.model = GDINA.options$higher.order.model,higher.order.method = GDINA.options$higher.order.method,
                      verbose = FALSE, catprob.parm = GDINA.obj$catprob.parm,
                      mono.constraint = GDINA.options$mono.constraint,
                      empirical = GDINA.options$empirical, att.prior = GDINA.options$att.prior,
                      att.str = GDINA.options$att.str,
                      nstarts = 3, conv.crit = GDINA.options$conv.crit,
                      maxitr = GDINA.options$maxitr, higher.order.struc.parm =  GDINA.options$higher.order.parm)
    out.list[[r]] <- list(itemprob = unlist(boot.out$catprob.parm),delta = unlist(boot.out$delta.parm), HO = c(boot.out$higher.order.struc.parm))
  }
}

  itemprob <- delta <- HO.slope <- HO.intercept <- NULL
  # index for item probability (IP) and delta (DT)
  indexIP <- GDINA.obj$catprob.parm
  indexDT <- GDINA.obj$delta.parm
  for (j in 1:J){
    indexIP[[j]] <- rep(j,length(indexIP[[j]]))
    indexDT[[j]] <- rep(j,length(indexDT[[j]]))
  }
  indexIP <- unlist(indexIP) #finalized
  indexDT <- unlist(indexDT) #finalized

  for (r in 1:Rsample){
    itemprob <- rbind(itemprob,out.list[[r]]$itemprob)
    delta <- rbind(delta,out.list[[r]]$delta)
    if (GDINA.options$higher.order) {
      HO.slope <- rbind(HO.slope,out.list[[r]]$HO[1:K])
      HO.intercept <- rbind(HO.intercept,out.list[[r]]$HO[(K+1):(2*K)])
    }

  }
  IP.var <- var(itemprob)
  DT.var <- var(delta)
  IP.SE <- sqrt(diag(IP.var))
  DT.SE <- sqrt(diag(DT.var))
  names(DT.SE) <- NULL

  IP.SE.list <- GDINA.obj$catprob.parm
  DT.SE.list <- GDINA.obj$delta.parm
  #SE in list format
  for(j in 1:J){
    IP.SE.list[[j]] <- IP.SE[which(indexIP==j)]
    DT.SE.list[[j]] <- DT.SE[which(indexDT==j)]
    names(IP.SE.list[[j]]) <- paste("P(",apply(attributepattern(log(length(which(indexIP==j)
    )) / log(2)),1,function(x){paste(x,collapse = "")}),")",sep = "")
  }
  intercept.var <- slope.var <- HO.intercept.resamples <-  NULL
    HO.slope.resamples <-  HO.slope.var <-  HO.intercept.var <-  NULL
    HO.intercept.se = HO.slope.se  <-  NULL
  if(GDINA.options$higher.order){
    slope.var <- var(HO.slope)
  intercept.var <- var(HO.intercept)
  HO.SE <- sqrt(data.frame(slope.SE=diag(slope.var),intercept.SE=diag(intercept.var)))
  HO.slope.se = round(HO.SE$slope.SE,digits)
  HO.intercept.se = round(HO.SE$intercept.SE,digits)
  HO.slope.var = round(slope.var,digits)
  HO.intercept.var = round(intercept.var,digits)
  HO.slope.resamples = round(HO.slope,digits)
  HO.intercept.resamples =round(HO.intercept,digits)
  }

  return(list(item.prob.se = lapply(IP.SE.list,round,digits),delta.se = lapply(DT.SE.list,round,digits),
              HO.slope.se = HO.slope.se,
              HO.intercept.se = HO.intercept.se ,
              item.prob.cov =round(IP.var,digits), delta.cov = round(DT.var,digits),HO.slope.var = slope.var,
              HO.intercept.var = HO.intercept.var,
              item.prob.resamples = round(itemprob,digits), delta.resamples = round(delta,digits),
              HO.slope.resamples = HO.slope,
              HO.intercept.resamples =HO.intercept))
}




