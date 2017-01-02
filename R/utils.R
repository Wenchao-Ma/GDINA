# see the example in ?integer
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

is.nonNegativeInteger <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol & x >= 0

is.positiveInteger <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol & x > 0


inputcheck <- function(dat, Q, model, sequential,higher.order,
                       higher.order.model, higher.order.method,
                       verbose, catprob.parm,mono.constraint,
                       empirical, att.prior, lower.p, upper.p,att.str,
                       nstarts, conv.crit, maxitr, higher.order.parm,
                       digits,diagnosis){
  # if (!max(dat,na.rm = TRUE)==1 | !min(dat,na.rm = TRUE)==0 ) stop("Data matrix can only contain 0, 1 or NA.",call. = FALSE)
  if(!is.logical(sequential)) stop("sequential must be logical.",call. = FALSE)
  if (!all(is.nonNegativeInteger(Q))) stop("Q matrix can only contain 0 and positive integers.",call. = FALSE)
  if (!is.matrix(dat) & !is.data.frame(dat)) stop("Data must be a matrix or data frame.",call. = FALSE)
  if (!is.matrix(Q) & !is.data.frame(Q)) stop("Q-matrix must be a matrix or data frame.",call. = FALSE)
  if (any(rowSums(Q)<1)) stop("Some items do not require any attributes.",call. = FALSE)
  if (any(colSums(Q)<1)) stop("Some attributes are not required by any items.",call. = FALSE)
  if (!sequential&&ncol(dat)!=nrow(Q))stop("The number of columns in data must be equal to the number of rows in Q-matrix.",call. = FALSE)
  if (!is.logical(higher.order)) stop("higher.order must be TRUE or FALSE.",call. = FALSE)
  if (higher.order){
    if (!higher.order.model %in% c("Rasch","2PL","1PL")) stop("higher.order.model must be Rasch, 1PL or 2PL.",call. = FALSE)
    if (!higher.order.method %in%c("BL","MMLE")) stop("higher.order.method must be BL or MMLE.",call. = FALSE)

  }
  if (!verbose%in%c(0,1,2)) stop("verbose must be 0, 1, or 2.",call. = FALSE)
  if (!is.logical(att.str)) stop("att.str must be TRUE or FALSE.",call. = FALSE)
  if (!all(sapply(mono.constraint,is.logical))) stop("mono.constraint must be TRUE or FALSE.",call. = FALSE)
  if (!length(mono.constraint)%in%c(1,ncol(dat))) stop("length of mono.constraint must be equal to 1 or the number of items.",call. = FALSE)
  if (!is.logical(empirical)) stop("empirical must be TRUE or FALSE.",call. = FALSE)

  if (!is.positiveInteger(nstarts)) {nstarts <- 1; warning("nstarts must be a positive integer.")}
if (!is.null(catprob.parm)){
  if (!is.list(catprob.parm)) stop("catprob.parm must be a list.",call. = FALSE)
  if (length(catprob.parm)!=nrow(Q)) stop("The length of catprob.parm is not correct.",call. = FALSE)

}
  if (att.str) {
    if (max(Q)>1) stop("Attribute structure cannot be specified if attributes are polytomous.",call. = FALSE)
    if(higher.order) stop("Higher-order structure is not allowed if att.str = TRUE.",call.=FALSE)
  }


  if(any(lower.p>=upper.p)) stop("lower.p must be less than upper.p.",call. = FALSE)
  if(any(upper.p<0)||any(upper.p>1)) stop("upper.p must range from 0 to 1.",call. = FALSE)

  # if(is.null(att.prior) & higher.order == FALSE & empirical == FALSE) stop("When higher-order and empirical are FALSE, att.prior cannot be NULL.",call. = FALSE)
}


inputcheck.sim <- function(N, Q, gs.parm=NULL, sequential,model = "GDINA", type = "random",
                           catprob.parm = NULL, delta.parm = NULL)
  {
  if (!is.nonNegativeInteger(N) ) stop("N must be negative integer.",call. = FALSE)
  if (!all(is.nonNegativeInteger(Q))) stop("Q matrix can only contain 0 and positive integers.",call. = FALSE)
  if (!is.matrix(Q) & !is.data.frame(Q)) stop("Q-matrix must be a matrix or data frame.",call. = FALSE)
  if(!tolower(type) %in% c("random","equal")) stop("type must be either random or equal.",call. = FALSE)
  if(is.null(gs.parm)&&is.null(catprob.parm)&&is.null(delta.parm)) stop("Item parameters must be specified.",call. = FALSE)
  if (sum(c(!is.null(catprob.parm),!is.null(delta.parm),!is.null(gs.parm)))>1) stop("Item parameters can only be specified using one of itemprob.param, delta.param or gs.param.",call. = FALSE)
  if(!is.null(catprob.parm)&&!is.list(catprob.parm)) stop("itemprob.parm must be NULL or a list.",call. = FALSE)
  if(!is.null(gs.parm)&&!is.data.frame(gs.parm)&&!is.matrix(gs.parm)) stop("gs.parm must be NULL, a matrix or data frame.",call. = FALSE)
  if(!is.null(delta.parm)&&!is.list(delta.parm)) stop("delta.parm must be NULL or a list.",call. = FALSE)
  }

model.transform <- function(model,J){
  if(length(model)!=1&&length(model)!=J) stop("model must be a scalar or a vector with the same length as the test.", call. = FALSE)
  M <- c("GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM")
  if (is.character(model))
  {
    model <- toupper(model)
    if (!all(model %in% c("GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM")))
    {
      return(warning(
        "The model for each item can only be \"GDINA\",\"DINA\",\"DINO\",\"ACDM\",\"LLM\",or \"RRUM\"."
      ))
    }
    model <- match(model, M) - 1

  } else if (is.numeric(model))
  {
    if (!all(model %in% 0:5))
    {
      return(warning(
        "Model can only be \"GDINA\",\"DINA\",\"DINO\",\"ACDM\",\"LLM\",or \"RRUM\"."
      ))
    }
  }

  if (length(model) == 1)
  {
    model <- model * rep(1, J)
  } else if (length(model) != J)
  {
    return(warning("the length of model must be 1 or the same as the test length."))
  }
  return(model)
}


RDINA <- function(Ks){
  R <- list()
  for (k in Ks){
    if (k>1){
      nconstr <- 2^k-2 # number of constraints or number of rows of restriction matrix
      nzeros <- rep(0,nconstr - 1)
      R[[k]] <- cbind(matrix(c(1, rep(c(nzeros,-1,1),nconstr - 1),nzeros,-1),nrow = nconstr),0)
    }
  }
  return(R)
}


RDINO <- function(Ks){
  R <- list()
  for (k in Ks){
    if (k>1){
      nconstr <- 2^k-2 # number of constraints or number of rows of restriction matrix
      nzeros <- rep(0,nconstr - 1)
      R[[k]] <- cbind(0,matrix(c(1, rep(c(nzeros,-1,1),nconstr - 1),nzeros,-1),nrow = nconstr))
    }
  }
  return(R)
}


RACDM <- function(Ks){
  R <- list()
  for (k in Ks){
    if (k==1) next
  alp <- alpha(k)
  R[[k]] <- matrix(0,nrow = 2^k-(k+1),ncol = 2^k)
  for(r in 1:nrow(R[[k]])){
    R[[k]][r,c(1,(r+k+1))] <- 1
    loc <- which(alp[r+k+1,]==1)
    negloc1 <- loc[1]+1
    matchvec <- rep(0,k)
    matchvec[loc[-1]] <- 1
    negloc2 <- which(apply(alp,1,function(x)all(x==matchvec)))
    R[[k]][r,c(negloc1,negloc2)] <- -1
  }
  }
  return(R)
}

logit <- function(p){
  return (log(p/(1-p)))
}

inv.logit <- function(x){
  return(exp(x)/(1+exp(x)))
}

# generating design matrix for GDINA model - used in M-step
designM_GDINA <- function(Kjj){
  Mj <- designM(alpha(Kjj),3) # start from the Mj matrix of A-CDM
  if (Kjj>1){
    for (l in 2:Kjj){
      comb <- combn(2:(Kjj+1), l)
      Mj <- cbind(Mj,apply(comb,2,function(x){apply(Mj[,x],1,prod)}))
    }
  }
  return(Mj)
}

# calculate delta paramters from item probability
calc_delta <- function(itmpar,model,Kj,digits=4){
  delta <- vector("list",nrow(itmpar))
  for (j in 1:nrow(itmpar)){
    if (model[j]==0){
      Mj <- designM_GDINA(Kj[j])
      delta[[j]] <- c(solve(t(Mj)%*%Mj)%*%t(Mj)%*%itmpar[j,1:2^Kj[j]])
      if(Kj[j]==1){
        names(delta[[j]]) <- c("d0","d1")
      }else{
        names(delta[[j]]) <- c("d0",paste("d",unlist(lapply(apply(alpha(Kj[j]),1,function(x)which(x==1))[-1],function(x) paste(x,collapse = ""))),sep = ""))
        }
      }else if (model[j]==1|model[j]==2){
      delta[[j]] <- c(itmpar[j,1],itmpar[j,2^Kj[j]]-itmpar[j,1])
      names(delta[[j]]) <- c("d0","d1")
    }else if (model[j]>2){
      Mj <- designM(alpha(Kj[j]),model[j])
      if (model[j]==3) {
        delta[[j]] <- c(solve(t(Mj)%*%Mj)%*%t(Mj)%*%itmpar[j,1:2^Kj[j]])
      }else if(model[j]==4){
        delta[[j]] <- c(solve(t(Mj)%*%Mj)%*%t(Mj)%*%qlogis(itmpar[j,1:2^Kj[j]]))
      }else{
        delta[[j]] <- c(solve(t(Mj)%*%Mj)%*%t(Mj)%*%log(itmpar[j,1:2^Kj[j]]))
      }
      names(delta[[j]]) <- paste("d",0:Kj[j],sep = "")
    }
    delta[[j]] <- round(delta[[j]],digits)
  }
  return(delta)
}

format_delta <- function(delta,model,Kj,item.names = NULL,digits=4){
  #delta <- vector("list",nrow(itmpar))
  for (j in 1:length(delta)){
    if (model[j]==0){
      if(Kj[j]==1){
        names(delta[[j]]) <- c("d0","d1")
      }else{
        names(delta[[j]]) <- c("d0",paste("d",unlist(lapply(apply(alpha(Kj[j]),1,function(x)which(x==1))[-1],function(x) paste(x,collapse = ""))),sep = ""))
      }
    }else if (model[j]==1|model[j]==2){
      names(delta[[j]]) <- c("d0","d1")
    }else if (model[j]>2){
      names(delta[[j]]) <- paste("d",0:Kj[j],sep = "")
    }
    delta[[j]] <- round(delta[[j]],digits)
  }
  if (is.null(item.names)) item.names <- paste("Item", 1:length(delta))
  names(delta) <- item.names
  return(delta)
}

gs2p <- function(Q,
                 gs,
                 model,
                 type,
                 mono.constraint,
                 item.names=NULL,
                 digits = 4) {
  J <- nrow(Q)
  K <- ncol(Q)
  Kj <-rowSums(Q>0)  # The number of attributes for each item

  pattern <- alpha(K, T, Q)
  itemprob.matrix <- matrix(NA, J, 2 ^ max(Kj))
  L <- nrow(pattern)  # the number of latent groups
  #if(length(model)==1) M <- rep(model,J)
  delta.param <- itemprob.param <- vector("list", J)
  #calculate delta parameters in list format
  for (j in 1:J) {
    if (model[j] == 1 | model[j] == 2) {
      delta.param[[j]] <- c(gs[j, 1], 1 - gs[j, 2] - gs[j, 1])
    } else if (model[j] == 3) {
      p0 <- gs[j, 1]
      p1 <- 1 - gs[j, 2]
      if (type == "equal") {
        d <- rep((p1 - p0) / Kj[j], Kj[j])
      } else if (type == "random") {
        sumd <- p1 - p0
        if (Kj[j] == 1) {
          d <- sumd
        } else{
          d <- rep(0, Kj[j])
          for (k in 1:(Kj[j] - 1)) {
            d[k] <- runif(1, 0, sumd)
            sumd <- sumd - d[k]

          }
          d[Kj[j]] <- p1 - p0 - sum(d)
        }
      }
      delta.param[[j]] <- c(p0, d)

    } else if (model[j] == 4) {
      p0 <- plogis(gs[j, 1])
      p1 <- plogis(1 - gs[j, 2])
      if (type == "equal") {
        d <- rep((p1 - p0) / Kj[j], Kj[j])
      } else if (type == "random") {
        sumd <- p1 - p0
        if (Kj[j] == 1) {
          d <- sumd
        } else{
          d <- rep(0, Kj[j])
          for (k in 1:(Kj[j] - 1)) {
            d[k] <- runif(1, 0, sumd)
            sumd <- sumd - d[k]

          }
          d[Kj[j]] <- p1 - p0 - sum(d)
        }
      }
      delta.param[[j]] <- c(p0, d)
    } else if (model[j] == 5) {
      p0 <- log(gs[j, 1])
      p1 <- log(1 - gs[j, 2])
      if (tolower(type) == "equal") {
        d <- rep((p1 - p0) / Kj[j], Kj[j])
      } else if (tolower(type) == "random") {
        sumd <- p1 - p0
        if (Kj[j] == 1) {
          d <- sumd
        } else{
          d <- rep(0, Kj[j])
          for (k in 1:(Kj[j] - 1)) {
            d[k] <- runif(1, 0, sumd)
            sumd <- sumd - d[k]
            #print(sumd)
          }
          d[Kj[j]] <- p1 - p0 - sum(d)
        }
      }
      delta.param[[j]] <- c(p0, d)
    } else if (model[j] == 0) {
      p0 <- gs[j, 1]
      p1 <- 1 - gs[j, 2]
      if(Kj[j]==1){
        ps <- c(p0,p1)
      }else{
        if (mono.constraint[j]) {
          preloc <- preloclist(Kj[j])
          ps <- c(p0, rep(0, 2 ^ Kj[j] - 2), p1)
          for (l in 2:(length(ps) - 1)) {
            ps[l] <- runif(1, max(ps[preloc[[l]]]), p1)

          }
        } else{
          ps <- c(p0,runif(2 ^ Kj[j] - 2, 0, 1),p1)
        }
      }
      delta.param[[j]] <- c(solve(designM_GDINA(Kj[j])) %*% ps)
    }
  }



  for (j in 1:J) {
    Mj <- designmatrix(Kj[j], model[j])

    if (model[j] <= 3) {
      itemprob.matrix[j, 1:nrow(Mj)] <-
        itemprob.param[[j]] <- round(c(Mj %*% delta.param[[j]]), digits)
    } else if (model[j] == 4) {
      itemprob.matrix[j, 1:nrow(Mj)] <-
        itemprob.param[[j]] <-
        round(qlogis(c(Mj %*% delta.param[[j]])), digits)
    } else if (model[j] == 5) {
      itemprob.matrix[j, 1:nrow(Mj)] <-
        itemprob.param[[j]] <- round(exp(c(Mj %*% delta.param[[j]])), digits)
    }



    # prob[[j]] <- item.param[j,1:length(tmp)] <- round(tmp,digits)
    names(itemprob.param[[j]]) <-
        paste0("P(", apply(alpha(Kj[j]), 1, paste0, collapse=""), ")")


  }
  delta.param <- format_delta(delta.param, model, Kj, digits = digits)
  return(
    list(
      delta.parm = delta.param,
      itemprob.matrix = itemprob.matrix,
      itemprob.parm = itemprob.param
    )
  )
}

# generate which latent class should have a low probability success if monotonicity constraints are conformed
preloclist <- function(K){
  patt <- t(alpha(K))
  apply(patt, 2, function(x) {
    loc <- which(colSums((1-x)*(patt|x))==0)
    loc[-length(loc)]
    })
}



#list to matrix transformation
# modified from http://stackoverflow.com/questions/3699405/how-to-cbind-or-rbind-different-lengths-vectors-without-repeating-the-elements-o
l2m <- function(l, len=max(sapply(l,length)))
{
  t(sapply(l, 'length<-', value=len))
}

m2l <- function(m,remove=NA){
  if(is.na(remove)){
    lapply(seq_len(nrow(m)), function(i) m[i,!is.na(m[i, ]) ])
  }else{
    lapply(seq_len(nrow(m)), function(i) m[i,m[i, ]!=remove ])
  }

}

DS.obj <- function(d,K,model,prob){
  if(model<=3){
    exp.p <- c(designmatrix(K,model)%*%d)
  }else if(model==4){
    exp.p <- plogis(c(designmatrix(K,model)%*%d))
  }else if(model==5){
    exp.p <- exp(c(designmatrix(K,model)%*%d))
  }
  sum(abs(exp.p-prob)) #minimize
}

DS.const <- function(d,K,model,prob){
  if(model<=3){
    exp.p <- c(designmatrix(K,model)%*%d)
  }else if(model==4){
    exp.p <- plogis(c(designmatrix(K,model)%*%d))
  }else if(model==5){
    exp.p <- exp(c(designmatrix(K,model)%*%d))
  }
  exp.p
}

DS <- function(prob,model){
  K <- log2(length(prob))
  if(any(prob>1)||any(prob<0)) stop("prob must not be less than 0 and greater than 1.",call. = FALSE)
  if (!is.positiveInteger(K)) stop("The length of prob is not correct.",call. = FALSE)
  model <- model.transform(model,J=1)
  if(model==1||model==2){
    d <- c(0.2,0.6)
  }else if(model>=3){
    d <- c(0.1,rep(0.8/K,K))
  }else if(model==4){
    d <- c(qlogis(0.1),rep((qlogis(0.9)-qlogis(0.1))/K,K))
  }else if(model==5){
    d <- c(log(0.1),rep((log(0.9)-log(0.1))/K,K))
  }

   DSoptim <- Rsolnp::solnp(d,DS.obj,ineqfun = DS.const,ineqLB = rep(0,2^K),control=list(trace=0),
                            ineqUB = rep(1,2^K),K=K,model=model,prob=prob)
  exp.p <- DS.const(DSoptim$par,K,model,prob)
  DS <- DSoptim$values[length(DSoptim$values)]
  conv <- DSoptim$convergence
  return(list(DS=DS,exp.p=exp.p,conv=conv))
}

vec_mat_match <- function(v,m,dim){
  #v row vector
  #m matrix
  which(apply(m,dim,identical,v))
}

alpha <- function(K,poly=F,Q=NULL){

  if (!poly||max(Q)==1){ # --calculate dichotomous alpha patterns
    if (K==1){
      alpha <- matrix(c(0,1),ncol=1)
    }else{
      alpha <- diag(K)
      for (l in 2:K){
        alpha <- rbind(alpha,t(apply(combn(K,l),2,function(x){apply(alpha[x,],2,sum)})))
      }
      alpha <- rbind(0,alpha)
    }
  }else{#polytomous Q -- calculate polytomous alpha patterns -- Q matrix is required
    alpha <- expand.grid(lapply(apply(Q,2,max),seq,from=0))
  }
  colnames(alpha) <- paste("A",1:K,sep = "")
  return(alpha)
}

# # of latent classes
no_LC <- function(Q){
  prod(apply(Q, 2, function(x)
  {
    max(length(unique(x)),2)
  }))
}

which.max.randomtie <- function(x,na.rm=TRUE){
  loc <- which(x==max(x,na.rm = na.rm))
  if(length(loc)>1){
    loc <- sample(loc,1)
  }
  return(loc)
}

which.min.randomtie <- function(x,na.rm=TRUE){
  loc <- which(x==min(x,na.rm=na.rm))
  if(length(loc)>1){
    loc <- sample(loc,1)
  }
  return(loc)
}

seq_coding <- function(dat,Q){
  out <- NULL
  x=table(Q[,1])
  for (j in 1:ncol(dat)){for(s in 1:x[j]){
    tmp <- dat[,j]
    misind <- which(tmp<s-1,arr.ind = TRUE)
    ind1 <- which(tmp>=s,arr.ind = TRUE)
    ind0 <- which(tmp==s-1,arr.ind = TRUE)
    tmp[ind1] <- 1
    tmp[ind0] <- 0
    if(length(misind)>0) tmp[misind] <- NA

    out <- cbind(out,tmp)
  }}
  return(out)
}



bdiag <- function(mlist,fill=0){
  loc <- sapply(mlist,dim)
  out <- matrix(fill,rowSums(loc)[1],rowSums(loc)[2])
  cr <- cc <- 1
  for(i in 1:length(mlist)){
    out[cr:(cr+nrow(mlist[[i]])-1),cc:(cc+ncol(mlist[[i]])-1)] <- mlist[[i]]
    cr <- cr + nrow(mlist[[i]])
    cc <- cc + ncol(mlist[[i]])
  }
  out
}


delta_se <- function(object,type){
  Q <- internalextract(object,"Q")
  Kj <- rowSums(Q>0)
  pj <- itemparm(object)
  if(internalextract(object,"sequential")){
    Qc <- internalextract(object,"Qc")
    dat <- seq_coding(internalextract(object,"dat"),Qc)
  }else{
    dat <- internalextract(object,"dat")
  }
m <- internalextract(object,"models_numeric")
  scof <- scorefun(mX=as.matrix(dat),
                   mlogPost=as.matrix(internalextract(object,"logposterior.i")),
                   itmpar=as.matrix(internalextract(object,"catprob.matrix")),
                   parloc=eta.loc(internalextract(object,"Q")),
                   model=m)

  scorep <- scof$score
  ind <- scof$index + 1
  score <- se <- vector("list",internalextract(object,"ncat"))
  c2 <- NULL
  N <- extract(object,"nobs")
  for (j in 1:internalextract(object,"ncat")){
    scorepj <- as.matrix(scorep[,ind[which(ind[,2]==j),1]])
    if (m[j]<4){
      if(m[j]==1||m[j]==2) {
        score[[j]] <- scorepj%*%designmatrix(1,m[j])
      }else{
        score[[j]] <- scorepj%*%designmatrix(Kj[j],m[j])
      }

    }else if (m[j]==4){
      pw <- scorepj*outer(rep(1,N),pj[[j]]*(1-pj[[j]]))
      score[[j]] <- pw%*%designmatrix(Kj[j],3)
    }else if(m[j]==5){
      pw <- scorepj*outer(rep(1,N),pj[[j]])
      score[[j]] <- pw%*%designmatrix(Kj[j],3)
    }
c2 <- c(c2,rep(j,ncol(score[[j]])))
score[[j]][is.na(score[[j]])] <- 0
score[[j]] <- score[[j]] * as.numeric(!is.na(dat[,j]))
  }

  if(type == 1){
    vars <- bdiag(lapply(score,function(x) solve(crossprod(x))))
  }else if(type == 2){
    vars <- solve(crossprod(do.call(cbind,score)))
  }
  c1 <- 1:nrow(vars)
  se.c <- sqrt(diag(vars))
  for(j in 1:length(se)){
    se[[j]] <- se.c[c1[which(c2==j)]]
  }
  return(list(cov=vars,se=se,ind=data.frame(item=c2,loc=c1)))
}

# Only correct when model GDINA DINA or DINO
# For ACDM, LLM and RRUM, delta method needs to be used to calculate
# variance of item probabilities from delta parameters
itemprob_se <- function(object,type){
  Q <- internalextract(object,"Q")
  m <- model.transform(extract(object,"models"),nrow(Q))
  pj <- l2m(internalextract(object,what = "catprob.parm"))
  Lj <- 2^rowSums(Q>0)
  for(j in which(m %in% c(1,2))){#DINA or DINO
    pj[j,2] <- pj[j,Lj[j]]
      if (Lj[j]>2) pj[j,3:ncol(pj)] <- -1
  }
  if(internalextract(object,"sequential")){
    Qc <- internalextract(object,"Qc")
    dat <- seq_coding(internalextract(object,"dat"),Qc)
  }else{
    dat <- internalextract(object,"dat")
  }
vars <-SE(as.matrix(dat), as.matrix(extract(object,"logposterior.i")),
              as.matrix(pj), eta.loc(Q), m, as.matrix(1 - is.na(dat)), type)
    std.err <- vars$se
    for (j in which((m) %in% c(1, 2))) {
      std.err[j, Lj[j]] <- std.err[j, 2]
      std.err[j, 1:(Lj[j] - 1)] <- std.err[j, 1]
    }
  # std.err[std.err<0] <- NA
  se <- m2l(std.err,remove = -1)

  covIndex <- vars$index+1
  covs <- vars$Var
  return(list(cov=covs,se=se,ind=data.frame(item=covIndex[,2],loc=covIndex[,1])))
}

scorefunc <- function(object,...){
  if(internalextract(object,"sequential")){
    Qc <- internalextract(object,"Qc")
    dat <- seq_coding(internalextract(object,"dat"),Qc)
  }else{
    dat <- internalextract(object,"dat")
  }
  scof <- scorefun(mX=dat,
                   mlogPost=internalextract(object,"logposterior.i"),
                   itmpar=internalextract(object,"catprob.matrix"),
                   parloc=eta.loc(internalextract(object,"Q")),
                   model=internalextract(object,"models_numeric"))
  index = scof$index + 1
  colnames(index) <- c("Column","Cat","Parm")
  list(score = scof$score, index = index)
}

Rmatrix.vec <- function(K){
  patt <- alpha(K)
  eta <- eta.loc(patt[-c(1,nrow(patt)),])
  Rv <- vector("list",nrow(eta))
  for (r in 1:nrow(eta)){
    for(lc in seq_len(max(eta[r,]))){
      loc <- which(eta[r,]==lc)
      tmp <- matrix(0,length(loc)-1,2^K)
      tmp[,loc[1]] <- 1
      tmp[cbind(seq_len(length(loc)-1),loc[-1])] <- -1
      Rv[[r]] <- rbind(Rv[[r]],tmp)
    }
  }
  return(Rv)
}


Rmatrix.att <- function(K){
  R <- vector("list",K)
  Lk <- 2^K
  if (K<=1) return(warning("K must be 2 or more!"))
  # gives which groups should be set to equal
  pattK <- alpha(K)
  pattK_1 <- alpha(K-1)
  for(a in 1:K){
    Rk <- matrix(0,Lk/2,Lk)
    for (l in 1:nrow(pattK_1)){
      loc <- which(apply(pattK[,-a,drop=FALSE],1,function(x){all(x==pattK_1[l,])}))
      Rk[l,loc[1]] <- 1
      Rk[l,loc[2]] <- -1
    }
    R[[a]] <- Rk
  }
  return(R)


}


valQrate <- function(trueQ,misQ,valQ){
  Qs <- data.frame(trueQ=c(as.matrix(trueQ)),misQ=c(as.matrix(misQ)),valQ=c(as.matrix(valQ)))
  CR <- data.frame(true2mis=
                     apply(matrix(c(0,0,
                                    1,1,
                                    0,1,
                                    1,0),ncol = 2,byrow = TRUE),1,function(x)rowMatch(Qs[,-3],x)$count),
                   mis2val=apply(matrix(c(0,0,0,
                                          1,1,1,
                                          0,1,0,
                                          1,0,1),ncol = 3,byrow = TRUE),1,function(x)rowMatch(Qs,x)$count),
                   row.names = c("000/00","111/11","010/01","101/10"))
  return(CR)
}
