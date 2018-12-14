#' Estimating multiple-strategy cognitive diagnosis models
#'
#' An (experimental) function for calibrating the multiple-strategy CDMs for dichotomous response data (Ma & Guo, 2019)
#'
#' @param dat A required binary item response matrix
#' @param msQ A multiple-strategy Q-matrix; the first column gives item numbers and the second column gives
#' the strategy number. See examples.
#' @param model CDM used; can be \code{"DINA"},\code{"DINO"},\code{"ACDM"},\code{"LLM"}, and \code{"RRUM"}, representing
#' the GMS-DINA, GMS-DINO, GMS-ACDM, GMS-LLM and GMS-RRUM in Ma & Guo (2019), respectively. It can also be \code{"rDINA"} and \code{"rDINO"},
#' representing restricted GMS-DINA and GMS-DINO models where delta_{jm1} are equal for all strategies. Note that only a single model can be used
#' for the whole test.
#' @param s strategy selection parameter. It is equal to 1 by default.
#' @param att.prior mixing proportion parameters.
#' @param delta delta parameters in list format.
#' @param control a list of control arguments
#'
#' @return an object of class \code{GMSCDM} with the following components:
#' \describe{
#' \item{IRF}{A matrix of success probabilities for each latent class on each item (IRF)}
#' \item{delta}{A list of delta parameters}
#' \item{attribute}{A list of estimated attribute profiles including EAP, MLE and MAP estimates.}
#' \item{testfit}{A list of test fit statistics including deviance, number of parameters, AIC and BIC}
#' \item{sIRF}{strategy-specific item response function}
#' \item{pjmc}{Probability of adopting each strategy on each item for each latent class}
#' \item{sprv}{Strategy pravelence}
#' }
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu}}
#'
#' @seealso \code{\link{GDINA}} for MS-DINA model and single strategy CDMs,
#' and \code{\link{DTM}} for diagnostic tree model for multiple strategies in polytomous response data
#'
#'@references
#'
#' Ma, W., & Guo, W. (2019). Cognitive Diagnosis Models for Multiple Strategies. \emph{British Journal of Mathematical and Statistical Psychology.}
#'
#' @examples
#'\dontrun{
#' ##################
#' #
#' # data simulation
#' #
#' ##################
#' set.seed(123)
#' msQ <- matrix(
#' c(1,1,0,1,
#' 1,2,1,0,
#' 2,1,1,0,
#' 3,1,0,1,
#' 4,1,1,1,
#' 5,1,1,1),6,4,byrow = T)
#' # J x L - 00,10,01,11
#' LC.prob <- matrix(c(
#' 0.2,0.7727,0.5889,0.8125,
#' 0.1,0.9,0.1,0.9,
#' 0.1,0.1,0.8,0.8,
#' 0.2,0.5,0.4,0.7,
#' 0.2,0.4,0.7,0.9),5,4,byrow=TRUE)
#' N <- 10000
#' att <- sample(1:4,N,replace=TRUE)
#' dat <- 1*(t(LC.prob[,att])>matrix(runif(N*5),N,5))
#'
#'
#' est <- GMSCDM(dat,msQ)
#' # item response function
#' est$IRF
#' # strategy specific IRF
#' est$sIRF
#'
#'
#' ################################
#' #
#' # Example 14 from GDINA function
#' #
#' ################################
#' Q <- matrix(c(1,1,1,1,0,
#' 1,2,0,1,1,
#' 2,1,1,0,0,
#' 3,1,0,1,0,
#' 4,1,0,0,1,
#' 5,1,1,0,0,
#' 5,2,0,0,1),ncol = 5,byrow = TRUE)
#' d <- list(
#'   item1=c(0.2,0.7),
#'   item2=c(0.1,0.6),
#'   item3=c(0.2,0.6),
#'   item4=c(0.2,0.7),
#'   item5=c(0.1,0.8))
#'
#'   set.seed(123)
#' sim <- simGDINA(N=1000,Q = Q, delta.parm = d,
#'                model = c("MSDINA","MSDINA","DINA",
#'                          "DINA","DINA","MSDINA","MSDINA"))
#'
#' # simulated data
#' dat <- extract(sim,what = "dat")
#' # estimation
#' # MSDINA need to be specified for each strategy
#' est <- GDINA(dat,Q,model = c("MSDINA","MSDINA","DINA",
#'                              "DINA","DINA","MSDINA","MSDINA"),
#'              control = list(conv.type = "neg2LL",conv.crit = .01))
#'
#' # Approximate the MS-DINA model using GMS DINA model
#' est2 <- GMSCDM(dat, Q, model = "rDINA", s = 10,
#'                control = list(conv.type = "neg2LL",conv.crit = .01))
#' }
#'
#'
#' @export
#'


GMSCDM <- function(dat, msQ, model = "ACDM", s = 1, att.prior=NULL,
                   delta = NULL,
                   control = list()){

  # allow to change lower.p on 10/28
  if (exists(".Random.seed", .GlobalEnv)) {
    oldseed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- oldseed)
  }
  else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv))
  }
  myControl <- list(conv.crit = 0.0001, maxitr = 2000,
                    conv.type=c("ip","mp"),seed=123, lower.p = 1e-4,
                    solver = "nloptr",algorithm = "NLOPT_LD_SLSQP")
  control <- utils::modifyList(myControl,control)
  set.seed(control$seed)

  item.no <- msQ[,1]
  str.no <- msQ[,2]
  Q <- msQ[,-c(1:2)]

  N <- nrow(dat)

  J <- ncol(dat)

  S <- nrow(Q)

  K <- ncol(Q)

  Kj <- rowSums(Q)
  patt <- attributepattern(K)
  L <- nrow(patt)  # The number of latent groups

  sm <- vector("logical",J)
  for(j in 1:J) sm[j] <- sum(item.no==j)>1

  if(is.null(att.prior)){
    att.prior <- matrix(rep(1/L, L),ncol = 1)
    logprior <- matrix(log(att.prior),nrow = L,ncol = 1)
  }


  Qjf <- unrestrQ(msQ)
  Qjf <- Qjf[which(Qjf[,2]==1),]
  Qj <- Qjf[,-c(1:2)]
  models <- c("DINA","DINO","rDINA","rDINO","ACDM","LLM","RRUM")
  model_numeric <- which(model==models) - 2
  d <- list()
  if(is.null(delta)){

    for(j in 1:J){
      ub <- runif(1,0.75,0.99)
      lb <- runif(1,0.01,0.25)
      db <- ub - lb
      jrows <- which(item.no==j)
      qj <- apply(Q[jrows,,drop=FALSE],2,max)
      Kjmax <- which(qj==1)
      if(model_numeric==-1){
        d[[j]] <- c(lb,rep(db,sum(item.no==j)))
      }else if(model_numeric==0){
        d[[j]] <- c(lb,rep(db,sum(item.no==j)))
      }else if(model_numeric==1){
        d[[j]] <- c(lb,db)
      }else if(model_numeric==2){
        d[[j]] <- c(lb,db)
      }else if(model_numeric==3){
        dj <- db/Kj[jrows]
        dd <- Q[jrows,,drop=FALSE]
        dd[dd==0] <- NA
        for(rr in 1:nrow(dd))
          dd[rr,] <- dd[rr,]*dj[rr]
        d[[j]] <- c(lb,apply(dd[,Kjmax,drop=FALSE],2,min,na.rm=TRUE))
      }else if(model_numeric==4){
        dj <- (qlogis(ub)-qlogis(lb))/Kj[jrows]
        dd <- Q[jrows,,drop=FALSE]
        dd[dd==0] <- NA
        for(rr in 1:nrow(dd))
          dd[rr,] <- dd[rr,]*dj[rr]
        d[[j]] <- c(qlogis(lb),apply(dd[,Kjmax,drop=FALSE],2,min,na.rm=TRUE))
      }else if(model_numeric==5){
        dj <- (log(ub)-log(lb))/Kj[jrows]
        dd <- Q[jrows,,drop=FALSE]
        dd[dd==0] <- NA
        for(rr in 1:nrow(dd))
          dd[rr,] <- dd[rr,]*dj[rr]
        d[[j]] <- c(log(lb),apply(dd[,Kjmax,drop=FALSE],2,min,na.rm=TRUE))
      }
    }
  }else{
    d <- delta
  }

  des <- list()
  i <- 1
  for(j in 1:J){
    jrows <- which(item.no==j)
    qj <- apply(Q[jrows,,drop=FALSE],2,max)
    Kjmax <- which(qj==1)
    des1 <- patt[,Kjmax,drop=FALSE]
    for(s in str.no[jrows]){
      if(model_numeric>2){
        red.qjs <- Q[which(item.no==j&str.no==s),which(qj==1)]
        des[[i]] <- cbind(1,sweep(des1,2,red.qjs,FUN = "*"))
      }else if(model_numeric==1){
        des[[i]] <- designmatrix(K,model_numeric)
        des[[i]][which(apply(patt[,which(Q[which(item.no==j&str.no==s),]==1),drop=FALSE],1,min)==1),2] <- 1
      }else if(model_numeric==2){
        des[[i]] <- designmatrix(K,model_numeric)
        des[[i]][which(apply(patt[,which(Q[which(item.no==j&str.no==s),]==1),drop=FALSE],1,max)==1),2] <- 1
      }else if(model_numeric==-1){# DINA - delta1 varies
        if(sm[j]){
          des[[i]] <- cbind(1,matrix(0,L,length(jrows)))
          des[[i]][which(apply(patt[,which(Q[which(item.no==j&str.no==s),]==1),drop=FALSE],1,min)==1),s+1] <- 1
        }else{
          des[[i]] <- designmatrix(K,1)
          des[[i]][which(apply(patt[,which(Q[which(item.no==j&str.no==s),]==1),drop=FALSE],1,min)==1),2] <- 1
        }

      }else if(model_numeric==0){#DINO - delta1 varies
        if(sm[j]){
          des[[i]] <- cbind(1,matrix(0,L,length(jrows)))
          des[[i]][which(apply(patt[,which(Q[which(item.no==j&str.no==s),]==1),drop=FALSE],1,max)==1),s+1] <- 1
        }else{
          des[[i]] <- designmatrix(K,2)
          des[[i]][which(apply(patt[,which(Q[which(item.no==j&str.no==s),]==1),drop=FALSE],1,max)==1),2] <- 1
        }
      }

      i <- i+1
    }

  }

  plc <- d2p(d = d,des = des,msQ = msQ,model = model_numeric,r=s)


  maxchg <- 10L
  itr <- 0L
  parm0 <- list(ip = c(plc), prior = c(exp(logprior)), neg2LL = 0, delt = unlist(d))

  dif.parm <- list(ip = 0,
                   prior = 0,
                   neg2LL = 0,
                   delt = 0)

  ############### START of while loop ##########
  while (maxchg > control$conv.crit && itr < control$maxitr)
  {
    d0 <- unlist(d)

    estep <- LikNR_LC(as.matrix(plc),
                      as.matrix(dat),
                      matrix(logprior,ncol = 1),
                      rep(1,N),
                      rep(1,N),
                      0)
    dev1 <- estep$LL*(-2)

    for(j in 1:J){
      if(control$solver=="auglag"){
        d[[j]] <- alabama::auglag(par = d[[j]],fn=objf,hin = inf,j=j,des=des,msQ=msQ,Nj=estep$N[j,],Rj = estep$R[j,],
                                  model=model_numeric,r=s,control.outer = list(trace=FALSE,method="nlminb",kkt2.check=FALSE,eps=1e-6))$par
      }else if(control$solver=="nloptr"){
        gf <- function(x,j,des,msQ,Nj,Rj,model,r) nloptr::nl.grad(x, objf, heps = 1e-7,des=des,msQ=msQ,Nj=estep$N[j,],Rj = estep$R[j,],
                                                                  model=model_numeric,r=r,j=j)
        ingf <- function(x,j,des,msQ,Nj,Rj,model,r) nloptr::nl.jacobian(x, inf3, heps = 1e-7,des=des,msQ=msQ,Nj=estep$N[j,],Rj = estep$R[j,],
                                                                        model=model_numeric,r=r,j=j)

        if(model=="ACDM"){
          lb <- rep(control$lower.p,length(d[[j]]))
          ub <- rep(Inf,length(d[[j]]))
        }else if(model=="LLM"){
          lb <- c(qlogis(control$lower.p),rep(0,length(d[[j]])-1))
          ub <- rep(Inf,length(d[[j]]))
        }else if(model=="RRUM"){
          lb <- c(log(control$lower.p),rep(0,length(d[[j]])-1))
          ub <- c(0,rep(-1*log(control$lower.p),length(d[[j]])-1))
        }else{
          lb <- rep(0,length(d[[j]]))
          ub <- rep(1,length(d[[j]]))
        }

        d[[j]] <- nloptr::nloptr(d[[j]],eval_f=objf,eval_grad_f = gf,eval_g_ineq = inf3,
                                 eval_jac_g_ineq = ingf,lb = lb, ub = ub,
                                 opts = list("algorithm"=control$algorithm,xtol_rel = 1e-4,print_level=0),
                                 des=des,msQ=msQ,Nj=estep$N[j,],Rj = estep$R[j,],
                                 model=model_numeric,r=s,j=j)$solution


      }

    }
    plc <- d2p(d = d,des = des,msQ = msQ,model=model_numeric,r=s)


    logprior <- c(estep$logprior)
    parm1 <- list(ip = c(plc), prior = c(exp(logprior)), neg2LL = -2*estep$LL, delt = unlist(d))
    dif.parm <- list(ip = max(abs(parm1$ip-parm0$ip),na.rm = TRUE),
                     prior = max(abs(parm1$prior-parm0$prior),na.rm = TRUE),
                     neg2LL = parm0$neg2LL-parm1$neg2LL,
                     delt = max(abs(parm1$delt-parm0$delt),na.rm = TRUE))

    parm0 <- parm1
    itr <- itr + 1
    maxchg <- 0
    if(any(tolower(control$conv.type)=="ip"))  maxchg <- max(maxchg,dif.parm$ip)
    if(any(tolower(control$conv.type)=="delta"))  maxchg <- max(maxchg,dif.parm$delt)
    if(any(tolower(control$conv.type)=="mp"))  maxchg <- max(maxchg,dif.parm$prior)
    if(any(tolower(control$conv.type)=="neg2ll"))  maxchg <- max(maxchg,abs(dif.parm$neg2LL))

   cat('\rIteration =',itr,' Max change =',formatC(maxchg,digits = 5, format = "f"),
        ' Deviance =',formatC(-2*estep$LL,digits = 2, format = "f"))

  }

  estep <- LikNR_LC(as.matrix(plc),
                    as.matrix(dat),
                    matrix(logprior,ncol = 1),
                    rep(1,N),
                    rep(1,N),
                    0)
  neg2LL <- -2 * estep$LL
  npar <- length(unlist(d)) + L - 1
  AIC <- 2 * npar + neg2LL
  BIC <- neg2LL + npar * log(N)


  sIRF <- pjmc <- sprv <- list()
  for(j in seq_len(J)){
    x <- sp(j,d[[j]],des,msQ,model=model_numeric,r=s)
    sIRF[[j]] <- x$sIRF #strategy specific IRF Pj(x=1|alpha_c,m)
    pjmc[[j]] <- x$pjmc #Pj(m|alpha_c)
    sprv[[j]] <- colSums(c(exp(estep$logprior))*x$pjmc) # strategy prevalence : Pj(m) = \sum_c Pj(m|alpha_c) p(alpha_c)
  }

  item.names <- paste("Item",seq_len(J))
  LC.names <- apply(patt,1,paste0,collapse="")

  sprv <- l2m(sprv)
  rownames(plc) <- rownames(sprv) <- item.names
  colnames(plc) <- paste0("P(",LC.names,")")
  colnames(sprv) <- paste("Strategy",seq_len(max(msQ[,2])))

  for(j in seq_len(J)){
    rownames(sIRF[[j]]) <- rownames(pjmc[[j]]) <- LC.names
    colnames(sIRF[[j]]) <- colnames(pjmc[[j]]) <- paste("Strategy",c(msQ[which(msQ[,1]==j),2]))
  }


  att <- list(EAP = 1*((exp(estep$logpost) %*% patt) > 0.5000),
              MLE = patt[apply(estep$loglik,1,which.max.randomtie),],
              MAP = patt[apply(estep$logpost,1,which.max.randomtie),])
  ret <- list(prob=plc,delta=d,att=att,testfit=list(deviance=neg2LL,npar=npar,AIC=AIC,BIC=BIC),sIRF=sIRF,pjmc=pjmc,sprv=sprv)

  class(ret) <- "GMSCDM"
  invisible(ret)
}

