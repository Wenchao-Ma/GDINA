#' Multiple-choice models
#'
#' This function estimates the multiple-choice DINA model (de la Torre, 2009).
#'
#' @param dat A required \eqn{N \times J} data matrix of N examinees to J items. Values must be 1, 2,... representing nominal categories.
#' Missing values are currently not allowed.
#' @param Qc A required category and attribute association matrix.
#'     The first column gives the item number, which must
#'     be numeric (i.e., 1,2,...) and match the number of column in the data.
#'     The second column indicates the coded category of each item. The number in the second column must match
#'     with the number in the data, but if a category is not coded, it should not be included in the Q-matrix.
#'     Entry 1 indicates that the attribute is measured by the category, and 0 otherwise.
#'     Note that the MC-DINA model assumes that the category with the largest number of 1s is the key and that the coded distractors should
#'     allow to assign examinees uniquely.
#' @param model \code{MCDINA} only currently. Other MC models may be incorporated.
#' @param key a numeric vector giving the key of each item. See \code{Examples}.
#'     \code{NULL} by default indicating the coded category requiring the largest number of 1s is the key.
#' @param group Group membership vector if a multiple group model is considered. It must only contain whole numbers and start from 1. The
#'     length must be equal to the number of rows of the data.
#' @param conv.crit The convergence criterion for max absolute change in \code{conv.type} for two consecutive iterations.
#' @param conv.type convergence criteria; Can be \code{pr} or \code{LL},
#'    indicating category response function, or -2 times log-likelihood,respectively.
#' @param maxitr The maximum iterations allowed.
#' @param SE logical; estimating standard error of item parameters? Default is \code{FALSE}.
#' @examples
#'\dontrun{
#'  # check the format of the data
#'  # Entry 0 is not allowed
#'  head(sim10MCDINA$simdat)
#'
#'  #---------------------------------
#'  # check the format of the Q-matrix
#'  #---------------------------------
#'  # Take item 1 as an example:
#'  # category 2 has a q-vector (1,0,0)
#'  # category 1 has a q-vector (0,1,0)
#'  # category 4 has a q-vector (1,1,0)
#'  # category 3 is not included in the Q-matrix because it is not coded
#'  # the order of the coded categories in the Q-matrix doesn't matter
#'
#'  sim10MCDINA$simQ
#'   #     Item coded cat A1 A2 A3
#'   #        1         2  1  0  0
#'   #        1         1  0  1  0
#'   #        1         4  1  1  0
#'   #...
#'  est <- MCmodel(sim10MCDINA$simdat,sim10MCDINA$simQ)
#'  est
#'  est$testfit
#'
#'  #--------------------------------------
#'  # Distractors involving more attributes
#'  #--------------------------------------
#'  # some distractors may involve attributes that are not invovled by the key option
#'  # this is not allowed by the "original" MC-DINA (de la Torre, 2009) but is allowed
#'  # in the current implementation
#'
#'  # Users need to specify the key for each item to appropriate handle such an issue
#'  # Note item 1 below: category 1 is the key (as indicated in the key argument below)
#'  # The distractor (category 4) involves an attribute not included by the key option
#'
#'  Qc <- matrix(c(1,	1,	1,	1,	0,
#'                 1,	2,	0,	1,	0,
#'                 1,	3,	1,	0,	0,
#'                 1,	4,	1,	0,	1,
#'                 2,	1,	1,	0,	0,
#'                 2,	3,	1,	1,	0,
#'                 2,	2,	1,	1,	1,
#'                 3,	4,	1,	1,	1,
#'                 3,	2,	1,	1,	0,
#'                 3,	3,	0,	1,	1,
#'                 4,	1,	0,	1,	1,
#'                 4,	2,	0,	0,	1,
#'                 5,	1,	1,	0,	0,
#'                 6,	3,	0,	1,	0,
#'                 7,	2,	0,	0,	1,
#'                 8,	4,	1,	0,	0,
#'                 9,	1,	0,	1,	0,
#'                 10, 4,	0,	0,	1),ncol = 5,byrow = TRUE)
#'
#'  est2 <- MCmodel(sim10MCDINA2$simdat,Qc, key = c(1,2,4,1,1,3,2,4,1,4))
#'  est2
#'  est2$prob.parm
#'  est2$testfit
#'  est2$attribute
#'
#'
#'  ###############################
#'  # a multiple group model
#'  ###############################
#' est <- MCmodel(sim10MCDINA$simdat, sim10MCDINA$simQ, group = rep(1:2, 1500))
#' # log posterior of different groups
#' est$lik$logprior
#'
#' }
#' @return an object of class \code{MCmodel} with the following components:
#' \describe{
#' \item{prob.parm}{A list of success probabilities for each reduced latent class on each item (IRF)}
#' \item{prob.se}{A list of standard errors of item parameters}
#' \item{attribute}{A list of estimated attribute profiles including EAP, MLE and MAP estimates.}
#' \item{testfit}{A list of test fit statistics including deviance, number of parameters, AIC and BIC}
#' \item{R}{expected # of individuals in each latent group choosing each option}
#' \item{lik}{posterior probability}
#' \item{itr}{Total # of iterations}
#' }
#'
#' @references
#'
#' De La Torre, J. (2009). A cognitive diagnosis model for cognitively based multiple-choice options. \emph{Applied Psychological Measurement, 33}, 163--183.
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive Diagnosis Modeling. \emph{Journal of Statistical Software, 93(14)}, 1-26.
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu}}
#'
#' @seealso \code{\link{GDINA}} for G-DINA model
#' @export
#'
#'
#'
MCmodel <- function(dat, Qc, model = "MCDINA", key = NULL,
                    group = NULL,conv.crit = .001, maxitr=2000, conv.type="pr", SE=FALSE){

  s1 <- Sys.time()
  mcm.call <- match.call()
  if (exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- oldseed)
  }else{
    on.exit(rm(".Random.seed", envir = .GlobalEnv))
  }
  if(any(dat==0))
    stop("data must only contain positive integers, e.g., 1,2,...",call. = FALSE)
  if(any(is.na(dat)))
    stop("Missing data are not allowed.",call. = FALSE)


  if(missing(dat)) missingMsg(dat)
  if(missing(Qc)) missingMsg(Qc)

  Q <- Qc[,-c(1:2)]

  N <- nrow(dat)

  J <- ncol(dat)


  if(!is.null(group)){
    g <- unique(group)
    no.g <- length(g)
    if(length(group)!=nrow(dat))
      stop("Length of group must be equal to the sample size in the data.", call. = FALSE)
  }else{
    no.g <- 1
    group <- rep(1,N)
  }


  unique.code <- apply(dat,2,unique)
  if(is.list(unique.code)){
    C <- sapply(unique.code,length) # total # of categories including 0
  }else if(is.matrix(unique.code)){
    C <- rep(nrow(unique.code),ncol(unique.code))
  }

  S0 <- sum(C)
  item.no <- rep(1:J,times = C - 1)
  item.no.0 <- rep(1:J,times = C)

  K <- ncol(Q)

  L <- 2^K
  init <- eta.mcm(Qc, model = model, no.options = C)

  param <- pl(init, model = model)
  eta <- init$eta

  logprior <- log(matrix(1/L, nrow=L, ncol = no.g))



  dif <- LL.2 <- 1
  itr <- 0
  while(itr < maxitr){ # E-M algorithm

    p0 <- param

    p <- pl2pm(param,eta)
    if(any(p<0)||any(p>1))
      stop("some item success probabilities are not between 0 and 1.",call. = FALSE)
    # print(p)
    # calculate likelihood and posterior
    likepost <- Lik_DTM_MG(as.matrix(p),as.matrix(dat-1),C-1,logprior,group)

    #number of expected examinees getting each score in each latent class - including 0 score
    R1 <- Rljs_DTM(likepost$logpost,dat-1,C-1) #S0 x L

    for(j in seq_len(J)){
      param[[j]] <- ColNormalize(aggregateCol(R1[which(item.no.0==j),],eta[j,]))
    }


    LL.1 <- -2*likepost$LL

    if(tolower(conv.type)=="LL")  {
      dif <- abs(LL.1-LL.2)
    }else if(tolower(conv.type)=="pr"){
      dif <- max(abs(unlist(param)-unlist(p0)))
    }
    itr <- itr + 1
    cat('\rIter =',itr,' Max. abs. change =',formatC(dif,digits = 5, format = "f"),
        ' Deviance  =',formatC(LL.1,digits = 3, format = "f"),'                                                                                 ')

    if(dif<conv.crit)
      break
    LL.2 <- LL.1

    logprior <- likepost$logprior

  }




  #---------update posterior
  p <- pl2pm(param,eta)
  likepost <- Lik_DTM_MG(as.matrix(p),as.matrix(dat-1),C-1,logprior,group)

  #number of expected examinees getting each score in each latent class - including 0 score
  R1 <- Rljs_DTM(likepost$logpost,dat-1,C-1) #S0 x L

  sco <- NULL
  for(j in 1:J){
    sco <- c(sco,score.mcm.j(Xj = dat[,j],
                             parloc.j = eta[j,],
                             catprob.j = param[[j]],
                             logpost = likepost$logpost))
  }
  se <- list()
  if(SE){
    SEs <- sqrt(diag(inverse_crossprod(do.call(cbind,sco))))

    st <- 1
    for(j in 1:J){
      se[[j]] <- rbind(matrix(SEs[st:(st+length(param[[j]])-ncol(param[[j]])-1)],ncol = ncol(param[[j]]),byrow = TRUE),NA)
      st <- st+length(param[[j]])-ncol(param[[j]])
    }
  }

  #-----------------Test Fit information----------------#

  neg2LL <- -2*likepost$LL

  npar.item <- 0
  for(j in seq_len(J))
    npar.item <- npar.item + length(param[[j]]) - ncol(param[[j]])

  npar <- L - 1 + npar.item

  AIC <- 2*npar + neg2LL

  BIC <- neg2LL + npar*log(N)

  test_fit <- list(neg2LL=neg2LL,npar=npar,AIC=AIC,BIC=BIC)
  #---------------Attribute estimation-----------------#
  patt <- attributepattern(K)
  att <- list(EAP=((exp(likepost$logpost)%*%patt)>0.5)*1,
              MAP=patt[apply(likepost$logpost,1,which.max.randomtie),],
              MLE=patt[apply(likepost$loglik,1,which.max.randomtie),])
  names(param) <- paste("Item",1:J)
  if (SE) names(se) <- paste("Item",1:J)
  total.row.label <- NULL
  for(j in seq_len(J)){
    colnames(param[[j]]) <- paste0("P(",init$label[[j]],")")
    rownames(param[[j]]) <- paste("Cat",seq_len(nrow(param[[j]])))
    if(SE){
      rownames(se[[j]]) <- paste("Cat",seq_len(nrow(param[[j]])))
      colnames(se[[j]]) <- paste0("SE(",init$label[[j]],")")
    }

    total.row.label <- c(total.row.label,paste("Item",j,"Cat",seq_len(nrow(param[[j]]))))
  }

  rownames(p) <- rownames(R1) <- total.row.label
  colnames(p) <- colnames(R1) <- colnames(eta)
  s2 <- Sys.time()

  ret <- list(prob.parm=param,prob.se=se,LC.prob = p,attribute=att,R=R1,lik=likepost,
              testfit=test_fit, time.used=s2-s1, call = mcm.call, dat=dat, Qc=Qc,
              itr=itr)

  class(ret) <- "MCmodel"

  return (ret)
}

eta.mcm <- function(Qc, model = "MCDINA",no.options = 4, key = NULL){
  Q <- Qc[,-c(1:2)]

  item.no <- Qc[,1]
  coded.cat.no <- Qc[,2]

  J <- length(unique(item.no))

  Cj <- tabulate(item.no)

  K <- ncol(Q)

  if(length(no.options)==1){
    no.options <- rep(no.options, J)
  }else if(length(no.options)!=J){
    stop("no.options must be of length 1 or the number of items.",call. = FALSE)
  }

  et <- matrix(0,J,2^K)
  et.label <- list()
  m <- list()
  if(model=="MCDINA"){
    for(j in 1:J){

      Qj <- Q[which(item.no==j),,drop=FALSE]
      Kj <- rowSums(Qj)
      kj.order <- order(Kj,decreasing = TRUE)
      coded.cat.no.j <- coded.cat.no[which(item.no==j)]
      if(any(duplicated(coded.cat.no.j))||any(Kj==0))
        stop("Q-matrix for item",j,"is not correctly specified.",call. = FALSE)

      if (Cj[j]>1){
        if(!is.null(key)){
          if(key[j]%in%coded.cat.no.j){
            key.loc <- which(coded.cat.no.j==key[j])
            if(key.loc!=kj.order[1]){
              kj.order <- c(key.loc,setdiff(kj.order,key.loc))
            }
          }else{
            stop("Option key for item",j,"is not given in the Q-matrix.",call. = FALSE)
          }
        }
        Qj <- Qj[kj.order,]
        et.label[[j]] <- c(apply(Qj,1,paste,collapse=""),paste(rep("*",ncol(Qj)),collapse = ""))
        m[[j]] <- matrix(0,no.options[j],length(et.label[[j]]))
        m[[j]][matrix(c(coded.cat.no.j[kj.order],seq_len(Cj[j])),nrow = Cj[j])] <- 1
        tmp <- eta(Qj)
        eta.j <- apply(rbind(0,(tmp==apply(tmp,1,max))*1),2,which.max)
        max.j <- max(eta.j)
        eta.j <- eta.j - 1
        eta.j[eta.j==0] <- max.j
        et[j,] <- eta.j

      }else{
        et.label[[j]] <- c(paste(Qj,collapse=""),paste(rep("*",ncol(Qj)),collapse = ""))
        eta.j <- eta(Qj)
        et[j,] <- 2 - (eta.j==max(eta.j))
        m[[j]] <- matrix(0,no.options[j],length(et.label[[j]]))
        m[[j]][matrix(c(coded.cat.no.j[kj.order],seq_len(Cj[j])),nrow = Cj[j])] <- 1
      }

    }
  }else if(model=="GNDM"){
    new.Q <- unrestrQ(Qc)[which(!duplicated(item.no)),-c(1:2)]
    et <- eta(new.Q)
    Kj <- rowSums(new.Q)
    for(j in 1:J){
      Qj <- Q[which(item.no==j),which(new.Q[j,]==1),drop=FALSE]
      Kjc <- rowSums(Qj)
      coded.cat.no.j <- coded.cat.no[which(item.no==j)]
      if(any(duplicated(coded.cat.no.j))||any(Kjc==0))
        stop("Q-matrix for item",j,"is not correctly specified.",call. = FALSE)
      att <- attributepattern(Kj[j])
      et.label[[j]] <- apply(att,1,paste,collapse="")
      m[[j]] <- matrix(0,no.options[j],length(et.label[[j]]))
      if(is.null(key)){
        loc.key <- which.max(Kjc)
      }else{
        if(!(key[j]%in%coded.cat.no.j)){
          stop("Option key for item",j,"is not given in the Q-matrix.",call. = FALSE)
        }else{
          loc.key <- which(coded.cat.no.j==key[j])
        }
      }
      tmp <- c(att%*%Qj[loc.key,])
      key.col <- which(tmp==max(tmp))
      m[[j]][coded.cat.no.j[loc.key],key.col] <- 1
      if(length(coded.cat.no.j)>1){
        for(jj in seq_len(length(coded.cat.no.j))){
          if(jj!=loc.key){
            tmp <- c(att%*%Qj[jj,])
            new.col <- setdiff(which(tmp==max(tmp)),key.col)
            m[[j]][coded.cat.no.j[jj],new.col] <- 1
            key.col <- c(new.col,key.col)
          }
        }
      }
    }
  }


  colnames(et) <- apply(attributepattern(K),1,paste,collapse="")
  rownames(et) <- names(et.label) <- paste("Item",seq_len(J))

  return(list(m=m,eta=et,label=et.label))

}

pl <- function(eta.obj,model = "MCDINA"){
  if(model=="MCDINA"){
    p <- eta.obj$m
    for(j in seq_len(length(p))){
      p[[j]][p[[j]]==1] <- .8
      p[[j]][p[[j]]==0] <- .2/(nrow(p[[j]])-1)
      p[[j]][,ncol(p[[j]])] <- 1/nrow(p[[j]])

    }
  }else if(model=="GNDM"){
    p <- eta.obj$m
    for(j in seq_len(length(p))){
      cs <- colSums(p[[j]])
      p[[j]][,which(cs==0)] <- 1/nrow(p[[j]])
      p[[j]][p[[j]]==1] <- .8
      p[[j]][p[[j]]==0] <- .2/(nrow(p[[j]])-1)
    }
  }
  p
}

pl2pm <- function(plist,eta){
  S0 <- sum(sapply(plist,nrow))
  L <- ncol(eta)
  J <- nrow(eta)
  p <- matrix(0,S0,L)

  for(l in seq_len(L)){
    loc <- 1
    for(j in seq_len(J)){
      pj <- plist[[j]][,eta[j,l]]
      p[loc:(loc+length(pj)-1),l] <- pj
      loc <- loc + length(pj)
    }
  }
  p
}


score.mcm.j <- function(Xj,
                        parloc.j,
                        catprob.j,
                        logpost){

  Xj[is.na(Xj)] <- -1 # remove NA from the data
  post <- exp(logpost) # posterior N x 2^K
  score.p <- vector("list",nrow(catprob.j)-1)
  for(r in 1:length(score.p)){
    score.p[[r]] <- aggregateCol(post,parloc.j)*
      (outer((Xj==r),catprob.j[r,],"/")-outer((Xj==nrow(catprob.j)),(catprob.j[nrow(catprob.j),]),"/"))
  }
  return(score.p)
}
