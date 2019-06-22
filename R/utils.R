# see the example in ?integer
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

is.nonNegativeInteger <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol & x >= 0

is.positiveInteger <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol & x > 0

getnames <- function(x){deparse(substitute(x))}

missingMsg <- function(x){
  stop(paste0('\"', x, '\" argument is missing.'), call.=FALSE)
}


inputcheck <- function(dat, Q, model, sequential,att.dist,latent.var,
                       verbose, catprob.parm,mono.constraint,
                       att.prior, lower.p, upper.p,att.str,
                       nstarts, conv.crit, maxitr){
  if(!is.logical(sequential)) stop("sequential must be logical.",call. = FALSE)
  if (!all(is.nonNegativeInteger(Q))) stop("Q matrix can only contain 0 and positive integers.",call. = FALSE)
  if (!is.matrix(dat) & !is.data.frame(dat)) stop("Data must be a matrix or data frame.",call. = FALSE)
  if(any(apply(dat,2,function(x) length(unique(x)))==1)) stop("Some items have only one response category and cannot be estimated.",call. = FALSE)
  if (!is.matrix(Q) & !is.data.frame(Q)) stop("Q-matrix must be a matrix or data frame.",call. = FALSE)
  if (any(rowSums(Q)<1)) stop("Some items do not require any attributes.",call. = FALSE)
  if (any(colSums(Q)<1)) stop("Some attributes are not required by any items.",call. = FALSE)
  if (ncol(dat)!=nrow(Q))stop("The number of columns in data does not match with the number of rows in Q-matrix.",call. = FALSE)
  if (!verbose%in%c(0,1,2)) stop("verbose must be 0, 1, or 2.",call. = FALSE)
  # if (!is.logical(att.str)) stop("att.str must be TRUE or FALSE.",call. = FALSE)
  if (!all(sapply(mono.constraint,is.logical))) stop("mono.constraint must be TRUE or FALSE.",call. = FALSE)
  if (!length(mono.constraint)%in%c(1,nrow(Q))) stop("Length of mono.constraint must be equal to 1 or the number of categories.",call. = FALSE)
  if (!is.positiveInteger(nstarts)) {nstarts <- 1; warning("nstarts must be a positive integer.")}
  if (!is.null(catprob.parm)){
    if (!is.list(catprob.parm)) stop("catprob.parm must be a list.",call. = FALSE)
    if (length(catprob.parm)!=nrow(Q)) stop("The length of catprob.parm is not correct.",call. = FALSE)

  }
  if(all(model==6)&&min(table(Q[,1]))==1) stop("You sure all models are MSDINA model?",call. = FALSE)
  if (max(dat, na.rm = TRUE) > 1 & !sequential)
    stop("Maximum response is greater than 1 - set sequential = TRUE to fit a sequential model.", call. = FALSE)
  if(tolower(latent.var)!="att"){
    if(!all(model%in%c(0,1,2))) stop("Only Bug DINA, DINO and G-DINA models are available.",call. = FALSE)
    if(any(mono.constraint))
      stop("Monotonic constraint is not allowed for the bug DINA, DINO and G-DINA models.",call. = FALSE)
  }
  if (!is.null(att.str)) {
    if (max(Q)>1) stop("Attribute structure cannot be specified if attributes are polytomous.",call. = FALSE)
    if(any(att.dist=="higher.order")) stop("Higher-order structure is not allowed if att.str = TRUE.",call.=FALSE)
    if(any(att.dist=="independent")) stop("Independent structure is not allowed if att.str = TRUE.",call.=FALSE)
    if(any(att.dist=="loglinear")) stop("Loglinear structure is not allowed if att.str = TRUE.",call.=FALSE)
    }


  if(any(lower.p>=upper.p)) stop("lower.p must be less than upper.p.",call. = FALSE)
  if(any(upper.p<0)||any(upper.p>1)) stop("upper.p must range from 0 to 1.",call. = FALSE)
}


inputcheck.sim <- function(N, Q, gs.parm=NULL, sequential, model = "GDINA", type = "random",
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
  if(sequential){
    if(max(rowSums(Q[,-c(1:2)]))==0)
      stop("Some rows of the Q-matrix contain only 0s.",call. = FALSE)
  }else{
    if(max(rowSums(Q))==0)
      stop("Some rows of the Q-matrix contain only 0s.",call. = FALSE)
  }
  }

model.transform <- function(model,J){
  if(length(model)==1){
    model <- rep(model, J)
  }else if(length(model)!=J){
    stop("model must be a scalar or a vector with the same length as the test.", call. = FALSE)
  }

  M <- c("LOGGDINA","LOGITGDINA","UDF","GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM","MSDINA")
  if (all(is.character(model)))
  {
    model <- toupper(model)
    if (!all(model %in% M))
     stop("The model for each item can only be \"UDF\",\"GDINA\",\"logitGDINA\",\"logGDINA\",\"DINA\",\"DINO\",\"ACDM\",\"LLM\", \"RRUM\", or \"MSDINA\".")

    model <- match(model, M) - 4 # GDINA is 0

  } else if (all(is.numeric(model)))
  {
    if (!all(model %in% -3:6))
      stop("The model for each item can only be \"UDF\",\"GDINA\",\"logitGDINA\",\"logGDINA\",\"DINA\",\"DINO\",\"ACDM\",\"LLM\", \"RRUM\", or \"MSDINA\".")
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
  alp <- alpha2(k)
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


# calculate delta paramters from item probability
calc_delta <- function(item.parm,DesignMatrices,linkfunc){
  if(is.matrix(item.parm)){
    item.parm <- m2l(item.parm)
  }
  delta <- vector("list",length(item.parm))
  for (j in 1:length(item.parm)){
    delta[[j]] <- c(Calc_Dj(item.parm[[j]], designMj = DesignMatrices[[j]], linkfunc = linkfunc[j]))
}
  return(delta)
}

format_delta <- function(delta,model,Kj,item.names = NULL,digits=4){
  if(!is.numeric(model))
    model <- match(model, c("lOGGDINA","LOGITGDINA","UDF", "GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM", "MSDINA")) - 4
  for (j in 1:length(delta)){
    if (model[j]%in%c(-3,-2,0)){
      if(Kj[j]==1){
        names(delta[[j]]) <- c("d0","d1")
      }else{
        name <- c("d0",paste("d",unlist(lapply(apply(alpha2(Kj[j]),1,function(x)which(x==1))[-1],function(x) paste(x,collapse = ""))),sep = ""))
        if(length(name)==length(delta[[j]]))
            names(delta[[j]]) <- name
      }
    }else if (model[j]%in%c(1,2,6)){
      names(delta[[j]]) <- c("d0","d1")
    }else if (model[j]%in%c(3:5)){
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
                 digits = 8){

  if(!is.numeric(model))
    model <- match(toupper(model), c("LOGGDINA","LOGITGDINA","UDF", "GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM", "MSDINA")) - 4
  J <- nrow(Q)
  K <- ncol(Q)
  Kj <-rowSums(Q>0)  # The number of attributes for each item

  pattern <- attributepattern(Q=Q)
  itemprob.matrix <- matrix(NA, J, 2 ^ max(Kj))
  L <- nrow(pattern)  # the number of latent groups
  #if(length(model)==1) M <- rep(model,J)
  delta.param <- itemprob.param <- vector("list", J)
  #calculate delta parameters in list format
  for (j in 1:J) {
    if (model[j]%in%c(1,2)) {
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
      d0 <- qlogis(gs[j, 1])
      dsum <- qlogis(1 - gs[j, 2])
      if (type == "equal") {
        d <- rep((dsum - d0) / Kj[j], Kj[j])
      } else if (type == "random") {
        sumd <- dsum - d0
        if (Kj[j] == 1) {
          d <- sumd
        } else{
          d <- rep(0, Kj[j])
          for (k in 1:(Kj[j] - 1)) {
            d[k] <- runif(1, 0, sumd)
            sumd <- sumd - d[k]

          }
          d[Kj[j]] <- dsum - d0 - sum(d)
        }
      }
      delta.param[[j]] <- c(d0, d)
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
      delta.param[[j]] <- c(solve(designmatrix(Kj[j],model[j])) %*% ps)
    }
  }



  for (j in 1:J) {
    Mj <- designmatrix(Kj[j], model[j])

    if (model[j] <= 3&&model[j]>=0) {
      itemprob.matrix[j, 1:nrow(Mj)] <-
        itemprob.param[[j]] <- round(c(Mj %*% delta.param[[j]]), digits)
    } else if (model[j] == 4) {
      itemprob.matrix[j, 1:nrow(Mj)] <-
        itemprob.param[[j]] <-
        round(plogis(c(Mj %*% delta.param[[j]])), digits)
    } else if (model[j] == 5) {
      itemprob.matrix[j, 1:nrow(Mj)] <-
        itemprob.param[[j]] <- round(exp(c(Mj %*% delta.param[[j]])), digits)
    }



    # prob[[j]] <- item.param[j,1:length(tmp)] <- round(tmp,digits)
    names(itemprob.param[[j]]) <-
        paste0("P(", apply(alpha2(Kj[j]), 1, paste0, collapse=""), ")")


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
  patt <- t(alpha2(K))
  apply(patt, 2, function(x) {
    loc <- which(colSums((1-x)*(patt|x))==0)
    loc[-length(loc)]
    })
}

gs2p.DTM <- function(Q,
                 gs,
                 model,
                 type,
                 mono.constraint,
                 linkfunc="logit",
                 item.names=NULL,
                 digits = 8) {
  J <- nrow(Q)
  K <- ncol(Q)
  Kj <-rowSums(Q>0)  # The number of attributes for each item

  pattern <- GDINA::attributepattern(K)
  itemprob.matrix <- matrix(NA, J, 2 ^ max(Kj))
  L <- nrow(pattern)  # the number of latent groups
  #if(length(model)==1) M <- rep(model,J)
  delta.param <- itemprob.param <- vector("list", J)
  #calculate delta parameters in list format
  for (j in 1:J) {
    des <- designmatrix(Kj[j],model = model[j])
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
      p0 <- qlogis(gs[j, 1])
      p1 <- qlogis(1 - gs[j, 2])
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
      tmp <- ps
      if(linkfunc=="logit")  tmp <- qlogis(ps)

      delta.param[[j]] <- c(solve(designmatrix(Kj[j])) %*% tmp)
    }
  }



  for (j in 1:J) {
    Mj <- designmatrix(Kj[j], model[j])

    if (model[j] <= 3 & model[j]>0) {
      itemprob.matrix[j, 1:nrow(Mj)] <-
        itemprob.param[[j]] <- round(c(Mj %*% delta.param[[j]]), digits)
    } else if (model[j] == 4) {
      itemprob.matrix[j, 1:nrow(Mj)] <-
        itemprob.param[[j]] <-
        round(qlogis(c(Mj %*% delta.param[[j]])), digits)
    } else if (model[j] == 5) {
      itemprob.matrix[j, 1:nrow(Mj)] <-
        itemprob.param[[j]] <- round(exp(c(Mj %*% delta.param[[j]])), digits)
    }else if(model[j]==0){
      tmp <- round(c(Mj %*% delta.param[[j]]), digits)
      if(linkfunc=="logit") itemprob.param[[j]] <- plogis(tmp) else itemprob.param[[j]] <- tmp
    }



    # prob[[j]] <- item.param[j,1:length(tmp)] <- round(tmp,digits)
    names(itemprob.param[[j]]) <-
      paste0("P(", apply(GDINA::attributepattern(Kj[j]), 1, paste0, collapse=""), ")")


  }
  delta.param <- format_delta(delta.param, model, Kj, digits = digits)
  return(
    list(
      delta.parm = delta.param,
      # itemprob.matrix = itemprob.matrix,
      itemprob.parm = itemprob.param
    )
  )
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
    exp.p <- c(designM(K,model)%*%d)
  }else if(model==4){
    exp.p <- plogis(c(designM(K,model)%*%d))
  }else if(model==5){
    exp.p <- exp(c(designM(K,model)%*%d))
  }
  sum(abs(exp.p-prob)) #minimize
}

DS.const <- function(d,K,model,prob){
  if(model<=3){
    exp.p <- c(designM(K,model)%*%d)
  }else if(model==4){
    exp.p <- plogis(c(designM(K,model)%*%d))
  }else if(model==5){
    exp.p <- exp(c(designM(K,model)%*%d))
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

seq_coding <- function(dat,Qc=NULL,Kj=NULL){
  out <- NULL
  if(is.null(Kj)){
    x <- table(Qc[,1])
  }else{
    x <- Kj
  }

  for (j in 1:ncol(dat)){
    for(s in 1:x[j]){
    tmp <- dat[,j]
    misind <- which(tmp<s-1,arr.ind = TRUE)
    ind1 <- which(tmp>=s,arr.ind = TRUE)
    ind0 <- which(tmp==s-1,arr.ind = TRUE)
    tmp[ind1] <- 1
    tmp[ind0] <- 0
    if(length(misind)>0) tmp[misind] <- NA

    out <- cbind(out,tmp)
    }
    }
  return(as.matrix(out))
}



bdiag <- function(mlist,fill=0){
  len <- length(mlist)
  for(r in len:1) if(is.null(mlist[[r]])) mlist[r] <- NULL
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






Rmatrix.vec <- function(K){
  patt <- attributepattern(K)
  Eta <- eta(patt[-c(1,nrow(patt)),,drop=FALSE])
  Rv <- vector("list",nrow(Eta))
  for (r in 1:nrow(Eta)){
    for(lc in seq_len(max(Eta[r,]))){
      loc <- which(Eta[r,]==lc)
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
  pattK <- attributepattern(K)
  pattK_1 <- attributepattern(K-1)
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



crossprod.na <- function(x, y, val=0) {
  crossprod(replace(x, is.na(x), val),
            replace(y, is.na(y), val)
  )
}


inverse_crossprod <- function(x) {
  if(!is.null(x))  MASS::ginv(crossprod(x))
}


designmatrix.bug <- function(Kj,model=1){
  if (!is.positiveInteger(Kj)) stop('Kj must be positive integer.',call. = FALSE)
  if (model==0){
    tmp <- designmatrix(Kj)
  }else if (model==1){ #DINA
    tmp <- matrix(1,2^Kj,2)
    tmp[2^Kj,2] <- 0
  }else if (model==2){ #DINO
    tmp <- matrix(c(rep(1,2^Kj+1),rep(0,2^Kj-1)),2^Kj,2)
  }
  return(tmp)
}

partial_order2 <- function(Kjj,AlphaPattern=NULL){
  if(is.null(AlphaPattern)){
    alp <- attributepattern(Kjj)
  }else{
    alp <- AlphaPattern
  }

  alp <- cbind(c(1:nrow(alp)),rowSums(alp),alp)
  out <- NULL

  for(k in 1:max(alp[,2])){
    for(i in 1:sum(alp[,2]==k-1)){
      alpk_1 <- alp[alp[,2]==k-1,,drop=FALSE]
      alpk <-  alp[alp[,2]==k,,drop=FALSE]
      out <- rbind(out,cbind(alpk[(apply(alpk,1,function(x){all(x-alpk_1[i,]>=0)})),1],alpk_1[i,1]))
    }
  }
  colnames(out) <- c("l","s")
  return(out)
}



dots <- function(name, value, ...) {

  # this function is copied from dots() from dots (Collado-Torres, 2014)
  # see https://github.com/lcolladotor/dots/blob/master/R/dots.R

  args <- list(...)

  if(!name %in% names(args)) {

    ## Default value

    return(value)

  } else {

    ## If the argument was defined in the ... part, return it

    return(args[[name]])

  }

}


RNjj <- function(fitGDINA) {
  # Date (version): 02/03/2018
  # Author: Miguel A. Sorrel
  # Input:
  ## fitGDINA: reduced model (DINA, DINO, or ACDM) fitted with the GDINA package

  # Comments
  ## Computes N and R for each latent group
  ## This function is required by the LM.CDM function of the GDINA package

  res <- list()

  dat <- extract(fitGDINA,"dat")
  q.matrix <- as.matrix(extract(fitGDINA,"Q"))
  J <- nrow(q.matrix)
  K <- ncol(q.matrix)
  L <- 2^K

  pattern=attributepattern(K)
  eta <- matrix(NA,2^K,J)
  for (l in 1:2^K) {
    for (j in 1:J) {
      Kj=sum(q.matrix[j,])
      Lj=which((q.matrix[j,]==1)!=0)
      Qj=attributepattern(Kj)
      tmp_pattern=pattern[l,Lj]
      eta[l,j]=which((apply((matrix(1,2^Kj,1)%*%tmp_pattern==Qj),1,prod)==1)!=0)
    }
  }
  post.i <- exp(fitGDINA$technicals$logposterior.i)

  for (jj in 1:J) {
    res[[jj]] <- matrix(data = 0,nrow = 2,ncol = length(unique(eta[,jj])))
    row.names(res[[jj]]) <- c("R.jj", "N.jj")

    for (gg in 1:length(unique(eta[,jj]))) {
      res[[jj]][1, gg] <- sum(post.i[, which(eta[,jj] == gg)] * dat[, jj])
      res[[jj]][2, gg] <- sum(post.i[, which(eta[,jj] == gg)])
    }

  }

  return(res)
}

itemprob_se_M <- function(object,type){

  dat <- extract(object,"dat")
  q.matrix <- as.matrix(extract(object,"Q"))
  m <- rep(3, times = nrow(q.matrix))
  pj <- l2m(extract(object,what = "catprob.parm"))
  Lj <- 2^rowSums(q.matrix>0)
  vars <-SE(as.matrix(dat), as.matrix(extract(object,"logposterior.i")),
            as.matrix(pj), eta(q.matrix), m, as.matrix(1 - is.na(dat)), type)
  std.err <- vars$se
  std.err[std.err<0] <- NA
  se <- m2l(std.err,remove = -1)

  covIndex <- vars$index+1
  covs <- vars$Var
  return(list(cov=covs,se=se,ind=data.frame(item=covIndex[,2],loc=covIndex[,1])))
}



##############
#  DTM
##############

# For item j, calculate the category response function (probabilities of getting score h) from the probabilities of pseduo items (nodes)
# > Tmatrix
# column 1: observed item responses
# column 2,...,M: the transformation matrix specifying how these observed responses are obtained from pseduo items (U1,U2,...,UM)
#       [,1] [,2] [,3] [,4]
# [1,]    0    0   NA   NA
# [2,]    1    1    0   NA
# [3,]    2    2   NA    0
# [4,]    3    1    1   NA
# [5,]    3    2   NA    1
# nodes.pr must be a list of length M-1.
# The first element is a matrix (L x 3) with 3 columns giving P(U1=0|alpha_l),P(U1=1|alpha_l),P(U1=2|alpha_l)
# The second and third elements each must be a matrix with 2 columns giving P(U2=0|alpha_l),P(U2=1|alpha_l) and P(U3=0|alpha_l),P(U3=1|alpha_l)
# Returned Pj - L x # of unique observed responses matrix
NodesP2ObsPj <- function(nodes.pr,Tmatrixj){
  catp <- NULL
  obsResp <- Tmatrixj[,1]
  uniqeObsResp <- unique(obsResp)
  Pj <- matrix(1,nrow = nrow(nodes.pr[[1]]),ncol = length(uniqeObsResp))
  for(x in 1:length(uniqeObsResp)){ # all possible observed responses
    Tx <- Tmatrixj[which(obsResp==uniqeObsResp[x]),,drop=FALSE]
    if(nrow(Tx)==1){ # only one way to achieve score x
      Txx <- Tx[,-1,drop=FALSE]
      Uloc <- which(!is.na(Txx))
      for(loc in Uloc){ # which nodes contribute to the prob
        Pj[,x] <- Pj[,x]*nodes.pr[[loc]][,1+Txx[1,loc]]
      }
    }else{ # more than one way to achieve score x
      Txx <- Tx[,-1,drop=FALSE]
      tmp <- matrix(1,nrow = nrow(Pj),ncol = nrow(Txx))
      for(w in 1:nrow(Txx)){
        Uloc <- which(!is.na(Txx[w,,drop=FALSE]))
        for(loc in Uloc){ # which nodes contribute to the prob
          tmp[,w] <- tmp[,w]*nodes.pr[[loc]][,1+Txx[w,loc]]
        }
      }
      Pj[,x] <- rowSums(tmp)
    }

  }
  return(Pj)
}

# Qcj <- Qc[which(item.no==j),,drop=FALSE]
# delta is ? x L matrix
# returned prj is nodes.pr in NodesP2ObsPj functions
NodesV2NodesPj <- function(Qcj,delta){

  pseudo.no <- Qcj[,2]
  prj <- vector("list",length(unique(pseudo.no)))
  for(i in unique(pseudo.no)){#pseudo items
    tmp <- cbind(1,t(exp(delta[which(pseudo.no==i),,drop=FALSE])))
    prj[[i]] <- tmp/rowSums(tmp)
  }
  return(prj)
}

# output.list - only applicable when type is tree: TRUE -> output is a list ; FALSE -> a matrix
v2p <- function(v,Qc,type="cumulative",linkfunc="identity",Tmatrix=NULL,output.list=FALSE){
  #Input: v S x L matrix
  #Output: p S0 x L matrix
  item.no <- Qc[,1]
  J <- length(unique(item.no))
  p <- NULL # will be a S0 x L matrix
  if(type=="cumulative"){
    if(linkfunc=="identity"){
      for(j in unique(item.no)) p <- rbind(p,-1*apply(rbind(1,v[which(item.no==j),,drop=FALSE],0),2,diff))
    }else if(linkfunc=="logit"){
      for(j in unique(item.no)) p <- rbind(p,-1*apply(rbind(1,plogis(v[which(item.no==j),,drop=FALSE]),0),2,diff))
    }
  }else if(type=="adjacent"){
    if(linkfunc=="identity"){
      for(j in unique(item.no)) {
        tmp <- apply(rbind(1,v[which(item.no==j),,drop=FALSE]),2,cumprod)/apply(rbind(1,1-v[which(item.no==j),,drop=FALSE]),2,cumprod)
        p <- rbind(p,tmp/matrix(colSums(tmp),nrow = nrow(tmp),ncol = ncol(tmp),byrow = TRUE))
      }
    }else if(linkfunc=="logit"){
      for(j in unique(item.no)) {
        tmp <- apply(rbind(1,exp(v[which(item.no==j),,drop=FALSE])),2,cumprod)
        p <- rbind(p,tmp/matrix(colSums(tmp),nrow = nrow(tmp),ncol = ncol(tmp),byrow = TRUE))
      }
    }
  }else if(type=="sequential"){
    if(linkfunc=="identity"){
      for(j in unique(item.no)) {
        tmp1 <- apply(rbind(1,v[which(item.no==j),,drop=FALSE]),2,cumprod)
        tmp2 <- rbind(1-v[which(item.no==j),,drop=FALSE],1)
        p <- rbind(p,tmp1*tmp2)
      }
    }else if(linkfunc=="logit"){
      for(j in unique(item.no)) {
        tmp1 <- apply(rbind(1,plogis(v[which(item.no==j),,drop=FALSE])),2,cumprod)
        tmp2 <- rbind(1-plogis(v[which(item.no==j),,drop=FALSE]),1)
        p <- rbind(p,tmp1*tmp2)
      }
    }
  }else if(type=="nominal"){
    if(linkfunc=="identity"){
      stop("Nominal model must be defined under logit link function.",call. = FALSE)
    }else if(linkfunc=="logit"){
      for(j in unique(item.no)) {
        tmp <- rbind(1,exp(v[which(item.no==j),,drop=FALSE]))
        p <- rbind(p,tmp/matrix(colSums(tmp),nrow = nrow(tmp),ncol = ncol(tmp),byrow = TRUE))
      }
    }
  }else if(type=="tree"){
    if(is.null(Tmatrix)) stop("Tmatrix must be specified for a tree model.",call. = FALSE)

    if(linkfunc=="logit"){
      pr <- list()
      for(j in unique(item.no)) {
        nodeprj <- NodesV2NodesPj(Qcj=Qc[which(item.no==j),,drop=FALSE],
                                  delta = v[which(item.no==j),,drop=FALSE])
        pr[[j]] <- NodesP2ObsPj(nodeprj,Tmatrix[[j]]) # P(X=h|alpha_l)
      }
      if(!output.list) p <- t(do.call(cbind,pr))
    }
  }
  return(p)
}


v2pj <- function(vj,type="cumulative",linkfunc="identity",Tmatrixj=NULL,Qcj=NULL){
  #vi must be a Sj x L matrix

  if(type=="cumulative"){
    if(linkfunc=="identity"){
      p <- -1*apply(rbind(1,vj,0),2,diff)
    }else if(linkfunc=="logit"){
      p <- -1*apply(rbind(1,plogis(vj),0),2,diff)
    }
  }else if(type=="adjacent"){
    if(linkfunc=="identity"){
      tmp <- apply(rbind(1,vj),2,cumprod)/apply(rbind(1,1-vj),2,cumprod)
    }else if(linkfunc=="logit"){
      tmp <- apply(rbind(1,plogis(vj)),2,cumprod)/apply(rbind(1,1-plogis(vj)),2,cumprod)
    }
    p <- tmp/matrix(colSums(tmp),nrow = nrow(tmp),ncol = ncol(tmp),byrow = TRUE)
  }else if(type=="sequential"){
    if(linkfunc=="identity"){
      p <- apply(rbind(1,vj),2,cumprod)*rbind(1-vj,1)
    }else if(linkfunc=="logit"){
      p <- apply(rbind(1,plogis(vj)),2,cumprod)*rbind(1-plogis(vj),1)
    }
  }else if(type=="nominal"){
    if(linkfunc=="identity"){
      stop("Nominal model must be defined under logit link function.",call. = FALSE)
    }else if(linkfunc=="logit"){
      tmp <- rbind(1,exp(vj))
      p <- tmp/matrix(colSums(tmp),nrow = nrow(tmp),ncol = ncol(tmp),byrow = TRUE)

    }
  }else if(type=="tree"){
    if(is.null(Tmatrixj)) stop("Tmatrix must be specified for a tree model.",call. = FALSE)

    if(linkfunc=="logit"){
      nodeprj <- NodesV2NodesPj(Qcj=Qcj, delta = vj)
      p <- t(NodesP2ObsPj(nodeprj,Tmatrixj)) # P(X=h|alpha_l)
    }
  }
  return(p) # Sj0 x L
}

null <- function(x) x
first.not.zero <- function(x) {
  if(all(x==0)) {
    y=NULL
  }else{
    y=x[which(x!=0)]
    y=y[1]
  }
  return(y)
}

parmtrans <- function(mparj,mIndj){
  parj <- apply(mparj*mIndj,2,first.not.zero)[which(colSums(mIndj)>0)]
  if(nrow(mIndj)>1) {
    parj <- c(mparj[1,1],diff(mparj[,1]),parj[-1])
  }else{
    parj <- unlist(parj)
  }
  return(parj)
}

invparmtrans <- function(vparj,mIndj){
  ncat <- nrow(mIndj)
  tmp <- rep(1,ncol(mIndj))# 2^K
  if(ncat>1){
    tmp[which(colSums(mIndj)>0)] <- vparj[-c(1:(ncat-1))]
    v <- tmp*t(mIndj)
    v[1,] <- cumsum(vparj[1:ncat])
  }else{
    tmp[which(colSums(mIndj)>0)] <- vparj
    v <- tmp*t(mIndj)
  }
  return(v)
}

model2numeric <- function(model,J=1){
  if(is.numeric(model)){
    if(J!=1&&length(model)!=J)
      model <- rep(model, J)
  }else{
  M <- c("lOGGDINA","LOGITGDINA","UDF", "GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM", "MSDINA")
  if(J!=1&&length(model)!=J)
    model <- rep(model, J)
  model <- match(toupper(model), M) - 4
  }
  model
}

model2character <- function(model,J=1){
  if(is.character(model)){
    if(J!=1&&length(model)!=J)
      model <- rep(model, J)
  }else{
    M <- c("LOGGDINA","LOGITGDINA","UDF", "GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM", "MSDINA")
    if(J!=1&&length(model)!=J)
      model <- rep(model, J)
    model <- M[model + 4]
  }
  model
}

linkf.numeric <- function(linkfunc, model.vector){


  J <- length(model.vector)
  # identitiy link -> 1
  # logit link -> 2
  # log link -> 3

  if(is.null(linkfunc)){

    linkfunc <- model2linkfunc(model.vector)


  }else if (length(linkfunc) == 1){
    linkfunc <- which(tolower(linkfunc)==c("identity","logit","log"))
    linkfunc <- rep(linkfunc, J)
  } else {
    if (length(linkfunc) != J) stop("Length of linkfunc is not correct.", call. = FALSE)
    tmp <- linkfunc
    linkfunc <- rep(1,J)
    linkfunc[which(tolower(tmp)=="logit")] <- 2
    linkfunc[which(tolower(tmp)=="log")] <- 3
  }
  linkfunc
}

model2rule.j <- function(model.j){

  if(is.character(model.j)){
    switch(toupper(model.j),
           LOGGDINA = 0,
           LOGITGDINA = 0,
           GDINA = 0,
           DINA = 1,
           DINO = 2,
           LLM = 3,
           ACDM = 3,
           RRUM = 3,
           MSDINA = 4,
           UDF = -1)
  }else if(is.numeric(model.j)){
    if(model.j%in%c(0,-3,-2)){
      cr <- 0
    }else if(model.j%in%c(4,5)){
      cr <- 3
    }else if(model.j==6){
      cr <- 4
    }else{
      cr <- model.j
    }

  }

  # 0 -> saturated model; 1 ->DINA; 2 ->DINO; 3 ->additive model; 4 ->MS-DINA; -1-> UDF
}
model2rule <- function(model.vector){
  sapply(model.vector,model2rule.j)
  # 0 -> saturated model; 1 ->DINA; 2 ->DINO; 3 ->additive model; 4 ->MS-DINA; -1-> UDF
}

model2linkfunc.j <- function(model.j){
  if(is.character(model.j)){
    lf <- switch(toupper(model.j),
           LOGGDINA = 3,
           LOGITGDINA = 2,
           GDINA = 1,
           DINA = 1,
           DINO = 1,
           LLM = 2,
           ACDM = 1,
           RRUM = 3,
           MSDINA = 1)
  }else if(is.numeric(model.j)){
    if(model.j%in%c(0:3,6)){
      lf <- 1
    }else if(model.j%in%c(4,-2)){
      lf <- 2
    }else if(model.j%in%c(5,-3)){
      lf <- 3
    }
  }
lf
  # 0 -> saturated model; 1 ->DINA; 2 ->DINO; 3 ->additive model; 4 ->MS-DINA; -1-> UDF
}
model2linkfunc <- function(model.vector){
  sapply(model.vector,model2linkfunc.j)
  # 0 -> saturated model; 1 ->DINA; 2 ->DINO; 3 ->additive model; 4 ->MS-DINA; -1-> UDF
}

initials <- function(Q,nstarts=1,DesignMatrices,latent.var="att",att.str=NULL){
  J <- length(DesignMatrices)
  ret <- list()
  if(tolower(latent.var)=="att"){
    if(nstarts==1){
      par <- list()
      for(j in seq_len(J)){
        g <- runif(1,0.05,0.25)
        p <- runif(1,0.75,0.95)
        d <- c(g,(p-g)/(ncol(DesignMatrices[[j]])-1))
        par[[j]] <- c(DesignMatrices[[j]] %*%d)
      }
      ret <- l2m(par)
    }else if(nstarts==3&is.null(att.str)){
      ret <- list(rep(NA,nstarts))
      for(i in seq_len(nstarts)){
        ret[[i]] <- gs2p(Q=Q,gs=matrix(runif(nrow(Q)*2,0.05,0.25),ncol = 2),
                         model = rep(i,nrow(Q)),type = "equal",
                         mono.constraint = rep(TRUE,nrow(Q)))$itemprob.matrix
      }

    }else{
      for(i in 1:nstarts){
        par <- list()
        for(j in seq_len(J)){
          g <- runif(1,0.01,0.45)
          p <- runif(1,0.65,0.99)
          dd <- runif(ncol(DesignMatrices[[j]])-1)
          d <- c(g,(p-g)*dd/sum(dd))
          par[[j]] <- c(DesignMatrices[[j]] %*% d)
        }
        ret[[i]] <- l2m(par)
      }

    }
  }else if(tolower(latent.var)=="bugs"){
    if(nstarts==1){
      ip <- gs2p(Q=Q,gs=matrix(runif(nrow(Q)*2,0.05,0.25),ncol = 2),
                 model = rep(3,nrow(Q)),type = "equal",
                 mono.constraint = TRUE)$itemprob.parm
      ip <- lapply(ip,function(x)1-x)
      if(is.list(ip)){
        ret <- l2m(ip)
      }else if(is.matrix(ip)){
        ret <- t(ip)
      }

    }else{
      par <- list(rep(NA,nstarts))
      for(i in 1:nstarts){
        ip <- gs2p(Q=Q,gs=matrix(runif(nrow(Q)*2,0.05,0.25),ncol = 2),
                   model = rep(3,nrow(Q)),type = "random",
                   mono.constraint = TRUE)$itemprob.parm
        ip <- lapply(ip,function(x)1-x)
        if(is.list(ip)){
          ret[[i]] <- l2m(ip)
        }else if(is.matrix(ip)){
          ret[[i]] <- t(ip)
        }

      }

    }
  }
  return (ret)
}



##############
# GMS CDMs
#
##############
which.max.randomtie <- function(x,na.rm=TRUE){
  loc <- which(x==max(x,na.rm = na.rm))
  if(length(loc)>1){
    loc <- sample(loc,1)
  }
  return(loc)
}
p2p <- function(v,r){
  sum(v^(r+1))/sum(v^r)
}
d2p <- function(d,des,msQ,model,r){
  L <- 2^(ncol(msQ)-2)
  item.no <- msQ[,1]
  str.no <- msQ[,2]
  J <- length(d)
  S <- length(des)
  plc <- matrix(0,J,L)
  for (j in 1:J){
    jrows <- which(item.no==j)
    pj <- matrix(0,L,length(jrows))
    for(s in str.no[jrows]){
      rowloc <- which(item.no==j&str.no==s)
      if(model<=3){
        pj[,s] <- c(des[[rowloc]]%*%d[[j]])
      }else if(model==4){
        pj[,s] <- plogis(c(des[[rowloc]]%*%d[[j]]))
      }else if(model==5){
        pj[,s] <- exp(c(des[[rowloc]]%*%d[[j]]))
      }

    }
    plc[j,] <- rowSums(pj^(r+1))/rowSums(pj^r)
  }
  plc
}
dj2pj <- function(j,dj,des,msQ,model,r){
  L <- 2^(ncol(msQ)-2)
  item.no <- msQ[,1]
  str.no <- msQ[,2]
  jrows <- which(item.no==j)
  pj <- matrix(0,L,length(jrows))
  for(s in str.no[jrows]){

    rowloc <- which(item.no==j&str.no==s)
    if(model<=3){
      pj[,s] <- c(des[[rowloc]]%*%dj)
    }else if(model==4){
      pj[,s] <- plogis(c(des[[rowloc]]%*%dj))
    }else if(model==5){
      pj[,s] <- exp(c(des[[rowloc]]%*%dj))
    }

  }
  rowSums(pj^(r+1))/rowSums(pj^r)

}
# calculate strategy prevalence for item j
sp <- function(j,dj,des,msQ,model,r){
  L <- 2^(ncol(msQ)-2)
  item.no <- msQ[,1]
  str.no <- msQ[,2]
  jrows <- which(item.no==j)
  pj <- matrix(0,L,length(jrows))
  for(s in str.no[jrows]){

    rowloc <- which(item.no==j&str.no==s)
    if(model<=3){
      pj[,s] <- c(des[[rowloc]]%*%dj)
    }else if(model==4){
      pj[,s] <- plogis(c(des[[rowloc]]%*%dj))
    }else if(model==5){
      pj[,s] <- exp(c(des[[rowloc]]%*%dj))
    }

  }

  list(sIRF=pj,pjmc=pj^(r)/rowSums(pj^r))
}

objf <- function(dj,j,des,msQ,Nj,Rj,model,r){
  pj <- dj2pj(j,dj,des,msQ,model,r)
  pj[pj<.Machine$double.eps] <- .Machine$double.eps
  pj[pj>1-.Machine$double.eps] <- 1-.Machine$double.eps
  -1*sum(Rj*log(pj)+(Nj-Rj)*log(1-pj))
}
inf <- function(dj,j,des,msQ,Nj,Rj,model,r){
  desj <- do.call(rbind,des[which(msQ[,1]==j)])

  if(model<=3){
    p <- c(desj%*%dj)
  }else if(model==4){
    p <- plogis(c(desj%*%dj))
  }else if(model==5){
    p <- exp(c(desj%*%dj))
  }
  c(p,1-p) - 1e-6
}
inf3 <- function(dj,j,des,msQ,Nj,Rj,model,r){
  desj <- do.call(rbind,des[which(msQ[,1]==j)])

  if(model<=3){
    p <- c(desj%*%dj)
  }else if(model==4){
    p <- plogis(c(desj%*%dj))
  }else if(model==5){
    p <- exp(c(desj%*%dj))
  }
  c(1e-6 - p, p - (1 - 1e-6))
}
