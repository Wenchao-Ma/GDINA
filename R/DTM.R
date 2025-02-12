#' Experimental function for diagnostic multiple-strategy CDMs
#'
#' This function estimates the diagnostic tree model (Ma, 2018) for polytomous responses with multiple strategies. It is an experimental function, and will be further optimized.
#'
#'
#' @param dat A required \eqn{N \times J} data matrix of N examinees to J items. Missing
#'    values are currently not allowed.
#' @param Qc A required \eqn{J \times K+2} category and attribute association matrix,
#'     where J represents the number of items or nonzero categories and K represents the
#'     number of attributes. Entry 1 indicates that the attribute is
#'     measured by the item, and 0 otherwise. The first column gives the item number, which must
#'     be numeric and match the number of column in the data. The second column indicates the category number.
#' @param delta initial item parameters
#' @param conv.crit The convergence criterion for max absolute change in item parameters.
#' @param conv.type convergence criteria; Can be \code{pr},\code{LL} and \code{delta},
#'    indicating category response function, log-likelihood and delta parameters,respectively.
#' @param maxitr The maximum iterations allowed.
#' @param Tmatrix The mapping matrix showing the relation between the OBSERVED responses (rows) and the PSEDUO items (columns);
#' The first column gives the observed responses.
#' @examples
#'\dontrun{
#' K=5
#' g=0.2
#' item.no <- rep(1:6,each=4)
#' # the first node has three response categories: 0, 1 and 2
#' node.no <- rep(c(1,1,2,3),6)
#' Q1 <- matrix(0,length(item.no),K)
#' Q2 <- cbind(7:(7+K-1),rep(1,K),diag(K))
#' for(j in 1:length(item.no)) {
#'   Q1[j,sample(1:K,sample(3,1))] <- 1
#' }
#' Qc <- rbind(cbind(item.no,node.no,Q1),Q2)
#' Tmatrix.set <- list(cbind(c(0,1,2,3,3),c(0,1,2,1,2),c(NA,0,NA,1,NA),c(NA,NA,0,NA,1)),
#' cbind(c(0,1,2,3,4),c(0,1,2,1,2),c(NA,0,NA,1,NA),c(NA,NA,0,NA,1)),
#' cbind(c(0,1),c(0,1)))
#' Tmatrix <- Tmatrix.set[c(1,1,1,1,1,1,rep(3,K))]
#' sim <- simDTM(N=2000,Qc=Qc,gs.parm=matrix(0.2,nrow(Qc),2),Tmatrix=Tmatrix)
#' est <- DTM(dat=sim$dat,Qc=Qc,Tmatrix = Tmatrix)
#' }
#' @references
#'
#' Ma, W. (2018). A Diagnostic Tree Model for Polytomous Responses with Multiple Strategies. \emph{British Journal of Mathematical and Statistical Psychology.}
#'
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#'
#' @seealso \code{\link{GDINA}} for MS-DINA model and single strategy CDMs,
#' and \code{\link{GMSCDM}} for generalized multiple strategies CDMs for dichotomous response data
#' @export
DTM <- function(dat, Qc, delta = NULL, Tmatrix = NULL, conv.crit = 0.001, conv.type = "pr",maxitr = 1000){

  s1 <- Sys.time()
  DTMcall <- match.call()
  if (exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- oldseed)
  }else{
    on.exit(rm(".Random.seed", envir = .GlobalEnv))
  }
  if(missing(dat)) missingMsg(dat)
  if(missing(Qc)) missingMsg(Qc)
  type <- "tree"
  eq.const <- FALSE
  linkfunc <- "logit"

  Q <- Qc[,-c(1:2)]

  N <- nrow(dat)

  J <- ncol(dat)

  S <- nrow(Q)

  K <- ncol(Q)

  L <- 2^K

  item.no <- Qc[,1]
  if(type=="tree"){
    C <- sapply(Tmatrix,function(x)length(unique(x[,1]))-1)
  }else{
    C <- table(item.no)
  }
  # print(C)
  item.no.0 <- rep(1:J,times=C+1)

  dsg <- GDINA::designmatrix(K)
  patt <- GDINA::attributepattern(K)
  # which delta parameters need to be estimated (1) or fixed to 0 (0)
  # S x L
  Iestpar <- dsg[apply(Q,1,function(x) GDINA::rowMatch(df=patt,vec=x)$row.no),]

  if(is.null(delta)){
    delta <- Iestpar
    d <- gs2p.DTM(Q,gs=matrix(runif(2*S,0.05,0.25),S,2),model = rep(0,S),type = "equal",mono.constraint = rep(T,S),linkfunc = "logit")$delta.parm
    for(s in 1:S) delta[s,which(Iestpar[s,]==1)] <- d[[s]]
  }

  logprior <- log(rep(1/L,L))

  dif <- LL.2 <- 1
  itr <- 0
  while(dif > conv.crit & itr < maxitr){ # E-M algorithm

    d_copy <- delta

    #--------E step
    #probability of getting each score for each latent class
    # v - S x L matrix
    v <- t(apply(delta,1,function(x)dsg%*%x))
    # print(round(delta,2))
    p <- v2p(v,Qc,type=type,linkfunc=linkfunc,Tmatrix = Tmatrix) #S0 x L
    if(any(p<0)) stop("some item success probabilities are less than 0; check your initial values.",call. = FALSE)
    # print(p)
    # calculate likelihood and posterior
    likepost <- Lik_DTM(as.matrix(p),as.matrix(dat),C,logprior)

    #number of expected examinees getting each score in each latent class - including 0 score
    R0 <- Rljs_DTM(likepost$logpost,dat,C) #S0 x L

    LL.1 <- -2*likepost$LL
    dif.LL <- abs(LL.1-LL.2)
    if(tolower(conv.type)=="ll")  {
      dif <- dif.LL

      if (dif < conv.crit) break
    }
    LL.2 <- LL.1
    # print(dif.LL)
    opts <- Mstep_DTM(delta=delta,R0=R0,Iestpar=Iestpar,item.no=item.no,item.no.0=item.no.0,
                  K = K,type=type,linkfunc=linkfunc,eq.const = eq.const,Qc=Qc,Tmatrix=Tmatrix)
    delta <- opts$delta
    # print(delta)
    # print(delta)
    if(conv.type=="delta"){
      dif <- max(abs(d_copy - delta),na.rm = TRUE)
    }else if(conv.type == "pr") {
      dif <- max(abs(d_copy - delta),na.rm = TRUE)
    }
    cat('\rIter =',itr,' Max. abs. change =',formatC(dif,digits = 5, format = "f"),
        ' Deviance  =',formatC(-2*likepost$LL,digits = 3, format = "f"),'                                                                                 ')

    itr <- itr + 1


    att <- exp(likepost$logprior)
    # att[which(att<1e-2)] <- 1e-2
    att <- att/sum(att)
    logprior <- log(att)

  }




  #---------update posterior
  # v - S x L matrix
  v <- t(apply(delta,1,function(x)dsg%*%x))
  p <- v2p(v,Qc,type=type,linkfunc=linkfunc,Tmatrix = Tmatrix) #S0 x L

  # calculate likelihood and posterior
  likepost <- Lik_DTM(p,dat,C,logprior)




  #-----------------Test Fit information----------------#

  neg2LL <- -2*likepost$LL

  npar <- L - 1

  npar <- npar + sum(Iestpar)

  AIC <- 2*npar + neg2LL

  # AICc <- AIC + 2*npar*(npar+1)/(nrow(dat)-npar-1)

  BIC <- neg2LL + npar*log(nrow(dat))

  test_fit <- list(neg2LL=neg2LL,npar=npar,AIC=AIC,BIC=BIC)
  #---------------Attribute estimation-----------------#
  att <- list(EAP=((exp(likepost$logpost)%*%patt)>0.5)*1,
              MAP=patt[apply(likepost$logpost,1,which.max),])

  d <- list()
  for(s in 1:S) d[[s]] <- delta[s,which(Iestpar[s,]==1)]


  s2 <- Sys.time()

  ret <- list(delta=delta,red_delta = d,attribute=att,prob=p,lik=likepost,
              testfit=test_fit,time.used=s2-s1, call = DTMcall)

  class(ret) <- "DTM"

  return (ret)
}
