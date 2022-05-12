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



#' Generate hierarchical attribute structures
#'
#' This function can be used to generate hierarchical attributes structures, and
#' to provide prior joint attribute distribution with hierarchical structures.
#'
#' @param hierarchy.list a list specifying the hierarchical structure between attributes. Each
#'     element in this list specifies a DIRECT prerequisite relation between two or more attributes.
#'     See \code{example} for more information.
#' @param K the number of attributes involved in the assessment
#' @param Q Q-matrix
#' @param att.prob How are the probabilities for latent classes simulated? It can be \code{"random"} or \code{"uniform"}.
#'
#' @return att.str reduced latent classes under the specified hierarchical structure
#' @return impossible.latentclass impossible latent classes under the specified hierarchical structure
#' @return att.prob probabilities for all latent classes; 0 for impossible latent classes
#'
#' @seealso \code{\link{GDINA}}, \code{\link{autoGDINA}}
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @export
#'
#' @examples
#' \dontrun{
#' #################
#' #
#' # Leighton et al. (2004, p.210)
#' #
#' ##################
#' # linear structure A1->A2->A3->A4->A5->A6
#' K <- 6
#' linear=list(c(1,2),c(2,3),c(3,4),c(4,5),c(5,6))
#' att.structure(linear,K)
#'
#' # convergent structure A1->A2->A3->A5->A6;A1->A2->A4->A5->A6
#' K <- 6
#' converg <- list(c(1,2),c(2,3),c(2,4),
#'                c(3,4,5), #this is how to show that either A3 or A4 is a prerequisite to A5
#'                c(5,6))
#'att.structure(converg,K)
#'
#' # convergent structure [the difference between this one and the previous one is that
#' #                       A3 and A4 are both needed in order to master A5]
#' K <- 6
#' converg2 <- list(c(1,2),c(2,3),c(2,4),
#'                c(3,5), #this is how to specify that both A3 and A4 are needed for A5
#'                c(4,5), #this is how to specify that both A3 and A4 are needed for A5
#'                c(5,6))
#'att.structure(converg2,K)
#'
#' # divergent structure A1->A2->A3;A1->A4->A5;A1->A4->A6
#' diverg <- list(c(1,2),
#'                c(2,3),
#'                c(1,4),
#'                c(4,5),
#'                c(4,6))
#'att.structure(diverg,K)
#'
#' # unstructured A1->A2;A1->A3;A1->A4;A1->A5;A1->A6
#' unstru <- list(c(1,2),c(1,3),c(1,4),c(1,5),c(1,6))
#' att.structure(unstru,K)
#'
#' ## See Example 4 and 5 in GDINA function
#'}
att.structure <- function(hierarchy.list=NULL,K,Q,att.prob="uniform"){
  patt <- attributepattern(K=K,Q=Q)
  if (!is.null(hierarchy.list)){

    impos.id <- lapply(hierarchy.list,function(x){
      k <- length(x)
      if(max(x)>K) stop("Maximum element of hierarchy list cannot be greater than K.",call. = FALSE)
      which(rowSums(patt[,x[1:(k-1)],drop=FALSE])<patt[,x[k]])

    })


    if (is.list(impos.id)) {
      impos.id <- unique(c(unlist(impos.id)))
    }else{
      impos.id <- unique(c(impos.id))
    }
    pos.id <- setdiff(seq(1,nrow(patt)),impos.id)
    red.patt <- patt[-impos.id,]
  }else{
    red.patt <- patt
    pos.id <- seq(1,nrow(patt))
    impos.id <- NULL
  }
  prior <- rep(0,nrow(patt))
  if(tolower(att.prob)=="uniform"){
    prior[pos.id] <- 1/length(pos.id)
  }else if(tolower(att.prob)=="random"){
    prior[pos.id] <- runif(length(pos.id))
    prior <- prior/sum(prior)
  }
  names(prior) <- apply(patt,1,function(x){paste(x,collapse = "")})

  if(!is.null(impos.id)){
    permissible.att.prob <- prior[-impos.id]
  }else{
    permissible.att.prob <- prior
  }
  return(list(att.str=red.patt,impermissible.latentclass=sort(impos.id),att.prob=prior,permissible.att.prob=permissible.att.prob))
}



#' Generate design matrix
#'
#' This function generates the design matrix for an item
#'
#' @param Kj Required except for the MS-DINA model; The number of attributes required for item j
#' @param model the model associated with the design matrix; It can be "GDINA","DINA","DINO", "ACDM","LLM", "RRUM", "MSDINA", "BUGDINO", and "SISM".
#'        The default is "GDINA". Note that models "LLM" and "RRUM" have the same design matrix as the "ACDM".
#' @param Qj the Q-matrix for item j; This is required for "MSDINA", and "SISM" models; The number of rows is equal to the number
#'        of strategies for "MSDINA", and the number of columns is equal to the number of attributes.
#' @param no.bugs the number of bugs (or misconceptions). Note that bugs must be given in the last no.bugs columns.
#' @return a design matrix (Mj). See de la Torre (2011) for details.
#' @references
#'
#' de la Torre, J. (2011). The generalized DINA model framework. \emph{Psychometrika, 76}, 179-199.
#'
#'
#' @examples
#' \dontrun{
#' designmatrix(Kj = 2, model = "GDINA")
#' designmatrix(Kj = 3, model = "DINA")
#' msQj <- matrix(c(1,0,0,1,
#'                  1,1,0,0),nrow=2,byrow=TRUE)
#' designmatrix(model = "MSDINA",Qj = msQj)
#' }
#'
#' @export
#'
#'
#'
designmatrix <- function(Kj = NULL, model = "GDINA", Qj = NULL,no.bugs = 0) {

  cr <- model2rule.j(model)
  if (cr == 4) {
    # MSDINA
    # Kj is not necessary

    if (is.null(Qj) || max(Qj) > 1){
      stop("Qj is not correctly specified for the MS-DINA model.", call. = F)
    }else if(nrow(Qj)==1){
      m <- 1
      Kj <- sum(Qj)
    }else{
      Qj <- as.matrix(Qj)
      if (any(colSums(Qj) == 0))
        Qj <- Qj[, -which(colSums(Qj) == 0)]
      Ks <- rowSums(Qj)
      Kj <- sum(apply(Qj, 2, max))
      patt <- attributepattern(Kj)

      D <- matrix(c(rep(1, nrow(patt)), colSums(Qj %*% t(patt) == Ks) > 0), ncol = 2)
    }

  }else if(cr==5){#BUG-DINO
    if (!is.null(Qj)){
      Kj <- sum(Qj)
    }

    D <- designM(Kj, 2)
    D[,2] <- 1- D[,2]
  }else if(cr==6){#SISM
    if (is.null(Qj) || max(Qj) > 1){
      stop("Qj is not correctly specified for the SISM model.", call. = F)
    }

    Kj <- sum(Qj)
    if(no.bugs==0){
      stop("no.bugs is not correctly specified for the SISM model.", call. = F)
    }


    n.bug.j <- sum(utils::tail(Qj,no.bugs))
    n.att.j <- Kj - n.bug.j
    if(n.att.j==0){
      D <- designmatrix.bug(n.bug.j,2)
    }else if(n.bug.j==0){
      D <- designM(n.att.j,1)
    }else{
      patt <- attributepattern(Kj)
      col.att <- col.bug <- rep(0,nrow(patt))
      col.att[rowMatch(patt[,seq(n.att.j),drop=FALSE],rep(1,n.att.j))$row.no] <- 1
      col.bug[rowMatch(patt[,-seq(n.att.j),drop=FALSE],rep(0,n.bug.j))$row.no] <- 1
      D <- as.matrix(data.frame(1,col.att,col.bug,col.att*col.bug))

    }



  }else if(cr!=-1){
    D <- designM(Kj, cr)
  }else{
    stop("design matrix cannot be calculated.",call. = FALSE)
  }
  return(D)

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
#' @param att.str attribute structure. See \code{GDINA} for details.
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

LC2LG <-  function(Q,sequential = FALSE,att.str=NULL){
  if(sequential) Q <- Q[,-c(1:2)]
  K <- ncol(Q)
  if(is.null(att.str)){ # no structure
    AlphaPattern <- as.matrix(att.structure(hierarchy.list = att.str,K = K,Q = Q,att.prob="uniform")$`att.str`)
    parloc <- eta(Q)  #J x L

  }else if(is.matrix(att.str)){
    AlphaPattern <- att.str
    parloc <- eta(Q, AlphaPattern)  #J x L

  }else{
    AlphaPattern <- as.matrix(att.structure(hierarchy.list = att.str,K = K,Q = Q,att.prob="uniform")$`att.str`)
    parloc <- eta(Q, AlphaPattern)  #J x L

  }

  colnames(parloc) <- apply(AlphaPattern,1,function(x)paste0(x,collapse = ""))
  rownames(parloc) <- rownames(Q)
  parloc
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

  if(!all(sapply(mlist,is.matrix))){
    for(m in seq_len(length(mlist))) mlist[[m]] <- as.matrix(mlist[[m]])
  }

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
#' Only applicable to saturated model ofr joint attribute distribution
#'
#' @param object an object of class GDINA
#' @param parm Either \code{delta} or \code{prob} indicating score function for delta parameters and
#' success probabily parameters
#' @return a list where elements give the score functions for each item or category
#' @export
#'
#' @examples
#'\dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' score(fit)
#' }
#'
#'
score <- function(object,parm="delta"){
  if(parm=="delta"){
    score_d(object)
  }else{
    if(!all(extract(object,"models_numeric")%in%c(0,1,2)))
      stop("Score functions of success probabilities are only available for the G-DINA, DINA and DINO models.",call. = FALSE)
    score_p(object)
  }
}



