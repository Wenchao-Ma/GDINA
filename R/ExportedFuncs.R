#' Generate all possible attribute patterns
#'
#' This function generates all possible attribute patterns. The Q-matrix needs to be specified for polytomous attributes.
#'
#' @param K number of attributes
#' @param Q Q-matrix; required when Q-matrix is polytomous
#'
#' @return A \eqn{2^K\times K} matrix consisting of attribute profiles for \eqn{2^K} latent classes
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @export
#' @examples
#' attributepattern(3)
#'
#' q <- matrix(scan(text = "0 1 2 1 0 1 1 2 0"),ncol = 3)
#' q
#' attributepattern(Q=q)
#'
#' q <- matrix(scan(text = "0 1 1 1 0 1 1 1 0"),ncol = 3)
#' q
#' attributepattern(K=ncol(q),Q=q)

attributepattern <-  function(K,Q){
  if(missing(K)){
    if(!missing(Q)){
      K <- ncol(Q)
    }else{
      stop("The number of attribute or Q-matrix is needed.",call. = FALSE)
    }
  }
  if (missing(Q)){
    alpha <- alpha2(K)
  }else if(max(Q)==1){
    alpha <- alpha2(K)
  }else{
    alpha <- alphap(apply(Q,2,max))
  }
  colnames(alpha) <- paste("A",1:K,sep = "")
  return(alpha)
}


#' Transformation between latent classes and latent groups
#'
#' This function gives the equivalent latent classes which have the same category success
#' probabilities for each category or item.
#'
#' @param Q A required \eqn{J \times K} binary Q-matrix. J represents test
#'    length and K represents the number of attributes of this test.
#'    Entry 1 at row j and column k represents the \eqn{k^{th}} attribute
#'    is measured by item \eqn{j}, and 0 means item \eqn{j} does not
#'    measure attribute \eqn{k}.
#' @param sequential logical; whether the Q-matrix is a Qc-matrix for sequential models?
#' @return An item or category by latent class matrix. In the G-DINA model,
#'    if item j measures \eqn{Kj} attributes, \eqn{2^K} latent classes can
#'    be combined into \eqn{2^{Kj}} latent groups. This matrix gives
#'    which latent group each of \eqn{2^K} latent classes belongs to for each item.
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @export
#' @examples
#' attributepattern(3)
#'
#' q <- matrix(scan(text = "0 1 0 1 0 1 1 1 0"),ncol = 3)
#' q
#' LC2LG(Q = q)
#'

LC2LG <-  function(Q,sequential = FALSE){
  if(sequential) Q <- Q[,-c(1:2)]
  out <- eta(Q)
  colnames(out) <- apply(attributepattern(ncol(Q)),1,function(x)paste0(x,collapse = ""))
  rownames(out) <- rownames(Q)
  out
}



#' Create a block diagonal matrix
#'
#' @param mlist a list of matrices
#' @param fill value to fill the non-diagnoal elements
#'
#' @return a block diagonal matrix
#' @export
#' @seealso \code{bdiag} in \pkg{Matrix}
#' @examples
#'
#' m1 <- bdiagMatrix(list(matrix(1:4, 2), diag(3)))
#' m2 <- bdiagMatrix(list(matrix(1:4, 2), diag(3)),fill = NA)
#'
bdiagMatrix <- function(mlist,fill=0){

  loc <- sapply(mlist,dim)
  out <- matrix(fill,rowSums(loc)[1],rowSums(loc)[2])
  cr <- cc <- 1
  for(i in 1:length(mlist)){
    out[cr:(cr+nrow(mlist[[i]])-1),cc:(cc+ncol(mlist[[i]])-1)] <- as.matrix(mlist[[i]])
    cr <- cr + nrow(mlist[[i]])
    cc <- cc + ncol(mlist[[i]])
  }
  out
}


#' Combine R Objects by Columns
#'
#' Combine a sequence of vector, matrix or data-frame arguments by columns. Vector is treated as a column matrix.
#'
#' @param ... vectors or matrices
#' @param fill a scalar used when these objects have different number of rows.
#'
#' @return a data frame
#' @export
#' @seealso \code{\link[base]{cbind}}
#' @examples
#' cjoint(2,c(1,2,3,4),matrix(1:6,2,3))
#' cjoint(v1 = 2, v2 = c(3,2), v3 = matrix(1:6,3,2),
#'        v4 = data.frame(c(3,4,5,6,7),rep("x",5)),fill = 99)
#'
cjoint <- function(...,fill = NA){
  arg <- list(...)
  if(length(arg) == 0) {
    return(NULL)
  }else{
    maxrow <- max(sapply(arg,function(x){nrow(as.matrix(x))}))
    r.arg <- lapply(arg,function(x){
      df <- data.frame(x)
      nrow.x <- nrow(df)
      if(nrow.x < maxrow) {
        rbind(df,matrix(fill,maxrow-nrow.x,ncol(df),dimnames = list(NULL,names(df))))
        }else{
          data.frame(x)
        }})
    out <- do.call(cbind,r.arg)
    names(out) <- NULL
    return(out)
  }

}


#' Count the frequency of a row vector in a data frame
#'
#' @param df a data frame or matrix
#' @param vec the vector for matching
#' @return count the number of vector vec in the data frame
#' @return row.no row numbers of the vector vec in the data frame
#' @export
#'
#' @examples
#'
#' df <- data.frame(V1=c(1L,2L),V2=LETTERS[1:3],V3=rep(1,12))
#' rowMatch(df,c(2,"B",1))
#'
#'
#'
rowMatch <- function(df,vec=NULL){
  logicalvec <- apply(df,1,paste0,collapse="")==paste0(vec,collapse = "")
  return(list(count=sum(logicalvec),row.no=which(logicalvec)))
}

#' Unique values in a vector
#'
#' @param vec a vector
#' @return sorted unique values
#' @export
#'
#' @examples
#'
#' vec <- c(4,2,3,5,4,4,4)
#' unique_only(vec)
#' # see the difference from unique
#' unique(vec)
#'
#' vec <- letters[1:5]
#' unique_only(vec)
#'
#'
#'
#' @seealso \link[base]{unique}
#'
unique_only <- function(vec){
  vec <- sort(vec)
  unique(vec)[which(table(vec)==1)]
}


#' Generate unrestricted Qc matrix from an restricted Qc matrix
#'
#' @param Qc an restricted Qc matrix
#' @return an unrestricted Qc matrix
#' @export
#'
#' @examples
#'
#' Qc <- sim21seqDINA$simQc
#' Qc
#' unrestrQ(Qc)
#'
#'
unrestrQ <- function(Qc){
  for (j in unique(Qc[,1])){ # item j
    if(sum(Qc[,1]==j)>1)  Qc[which(Qc[,1]==j),2+which(colSums(Qc[which(Qc[,1]==j),-c(1:2)])>0)] <- 1
  }
  return(Qc)
}


#' Score function
#'
#' Calculate score function for each dichotomous item or each nonzero category for polytomous items
#'
#' @param object an object of class GDINA
#' @param parm Either \code{delta} or \code{prob} indicating score function for delta parameters and
#' success probabily parameters
#' @return a list where elements give the score functions for each item or category
#' @export
#'
#'
#'
#'
score <- function(object,parm="delta"){
  if(parm=="delta"){
    score_d(object)
  }else{
    if(any(extract(object,"models_numeric")>2)) stop("Score functions of success probabilities cannot be calculated for ACDM, LLM and RRUM.",call. = FALSE)
    score_p(object)
  }
}

#' Functions for internal use
#'
#' @param ... arguments passed to internal functions
#' @name internal_GDINA
NULL


#' @name internal_Lik
#' @export
#' @rdname internal_GDINA
internal_Lik <- function(...){
args <- list(...)
stopifnot(length(args)>0)
if(is.null(args$group)) args$group <- rep(1,nrow(args$dat))
if(is.null(args$weights)) args$weights <- rep(1,nrow(args$dat))
if(is.null(args$simplify)) args$simplify <- 0
LikNR(mpar = args$catprob.parm,
      mX = args$dat,
      vlogPrior = args$logprior,
      vgroup = args$group,
      mloc = args$eta,
      weights = args$weights,
      simplify = args$simplify)
}

#' @name internal_aggregateCol
#' @export
#' @rdname internal_GDINA
internal_aggregateCol <- function(...){
  args <- list(...)
  stopifnot(length(args)>0)
  aggregateCol(mX=args$x,
               ind = args$by)
}

#' @name internal_eta
#' @export
#' @rdname internal_GDINA
internal_eta <- function(...){
  args <- list(...)
  stopifnot(length(args)>0)
  eta(Q=args$Q)
}

#' @name internal_matchMatrix
#' @export
#' @rdname internal_GDINA
internal_matchMatrix <- function(...){
  args <- list(...)
  stopifnot(length(args)>0)
  matchMatrix(A=args$A,B=args$B)
}

#' @name internal_RowNormalize
#' @export
#' @rdname internal_GDINA
internal_RowNormalize <- function(...){
  args <- list(...)
  stopifnot(length(args)>0)
  RowNormalize(args$x)
}

#' @name internal_ColNormalize
#' @export
#' @rdname internal_GDINA
internal_ColNormalize <- function(...){
  args <- list(...)
  stopifnot(length(args)>0)
  ColNormalize(args$x)
}

#' @name internal_RowProd
#' @export
#' @rdname internal_GDINA
internal_RowProd <- function(...){
  args <- list(...)
  stopifnot(length(args)>0)
  rowProd(args$x,args$v)
}

#' @name internal_l2m
#' @export
#' @rdname internal_GDINA
internal_l2m <- function(...){
  l2m(...)
}

#' @name internal_l2m
#' @export
#' @rdname internal_GDINA
internal_m2l <- function(...){
  m2l(...)
}

#' @name internal_uP
#' @export
#' @rdname internal_GDINA
internal_uP <- function(...){
  args <- list(...)
  stopifnot(length(args)>0)
  uP(mloc = args$eta, mpar = args$catprob.parm)
}
